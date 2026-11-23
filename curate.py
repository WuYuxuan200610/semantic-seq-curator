"""curate.py —— 把生物序列清洗成可直接进语言模型输入口的数据。

目标只有一个：**输出里不能有任何模型喂不进去的乱码**。
因此本工具是"目标字母表驱动"的：先声明模型输入口允许哪些字符，
再把每个非标准字符映射成该字母表里的合法 token（或删除），
最后断言输出 100% 合法；同时记录每个字符的语义与映射理由，便于审计。

字母表（可用 --alphabet 指定；默认 auto 按每行推断的类型选择）：
    dna              ACGT + N                      （DNA 语言模型常用词表）
    rna              ACGU + N
    nucleic-mixed    ACGTU + N
    dna-iupac        ACGT + 全部 IUPAC 模糊码 + N  （tokenizer 支持模糊码时用）
    rna-iupac        ACGU + 全部 IUPAC 模糊码 + N
    protein          20 种标准氨基酸 + X           （严格蛋白词表）
    protein-extended 20 + X + B Z J U O *          （ESM 类词表）

映射规则（保证输出字符必属目标字母表）：
    小写字母        → 大写（低置信度/可变碱基来源，计入审计）
    IUPAC 模糊码    → fallback（核酸 N / 蛋白 X）  ← 目标词表不支持时
    gap  '-' '.'    → fallback（默认）或删除（--gaps drop）
    占位符 *[k]*    → k 个 fallback（k 已实证为丢失碱基数）
    其它未知字符    → fallback（默认）或删除（--unknown drop）

说明：严格字母表下把 IUPAC 模糊码映射为 N，与 seismic cleanfa 的选择一致；
本工具的差别在于**显式声明目标字母表、保证 100% 合法、并逐字符记录
"原来是什么 → 映射成什么 → 为什么"**，让映射可审计、可随 tokenizer 切换。
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import os
import re
import sys

PLACEHOLDER_RE = re.compile(r"\*\[ *(\d+)\]\*")   # *[k]*：k = 上游丢失的碱基个数
NUCLEIC_LETTERS = set("ACGTUNRYSWKMBDHVX")
DNA_BASES = set("ACGT")
RNA_BASES = set("ACGU")
NUCLEIC_BASES = set("ACGTU")
IUPAC_AMBIG = set("NRYWKSMBDHVX")
PROTEIN_BASES = set("ACDEFGHIKLMNPQRSTVWY")
PROTEIN_AMBIG = set("BJZX")
PROTEIN_MODIFICATIONS = set("UO*")
GAPS = set("-.")

TYPE_CHOICES = ("auto", "dna", "rna", "nucleic", "nucleic_mixed", "protein")
PLACEHOLDER_MODES = ("expand", "keep")
GAP_MODES = ("substitute", "drop")
UNKNOWN_MODES = ("substitute", "drop")

# 目标字母表：allowed = 可直接进模型的字符；fallback = 映射目标；kind 决定回退字符
ALPHABET_PROFILES = {
    "dna": {"allowed": set("ACGTN"), "fallback": "N", "kind": "nucleic"},
    "rna": {"allowed": set("ACGUN"), "fallback": "N", "kind": "nucleic"},
    "nucleic-mixed": {"allowed": set("ACGTUN"), "fallback": "N", "kind": "nucleic"},
    "dna-iupac": {"allowed": set("ACGTNRYWSKMBDHV"), "fallback": "N", "kind": "nucleic"},
    "rna-iupac": {"allowed": set("ACGUNRYWSKMBDHV"), "fallback": "N", "kind": "nucleic"},
    "protein": {"allowed": set("ACDEFGHIKLMNPQRSTVWYX"), "fallback": "X", "kind": "protein"},
    "protein-extended": {"allowed": set("ACDEFGHIKLMNPQRSTVWYXBZJUO*"),
                         "fallback": "X", "kind": "protein"},
}
ALPHABET_CHOICES = ("auto",) + tuple(ALPHABET_PROFILES)
DECISION_REASONS = {
    "小写→大写": "比对惯例中大小写代表低置信度/可变位点；大写是该字母的合法 token。",
    "模糊码→fallback": "目标字母表不含该模糊码；映射为 N/X（未知）以成为合法输入。",
    "模糊码(合法·保留)": "目标字母表（*-iupac）本身支持该模糊码，无需映射。",
    "gap→fallback": "目标字母表不含 gap 字符；映射为 N/X 或按 --gaps drop 删除。",
    "gap(删除)": "按 --gaps drop 删除，序列变为去比对形式（长度会变）。",
    "占位符→fallback×k": "上游丢失 k 个碱基（k 已 FASTA 侧翼锚定实证），还原为 k 个 N/X。",
    "未知→fallback": "目标字母表不含该字符；映射为 N/X 以保证可输入，不编造语义。",
    "未知(删除)": "按 --unknown drop 删除无法归类且不在字母表内的字符。",
}


class CurateInputError(ValueError):
    """输入不足以安全处理；必须让用户知道，而不是静默猜测。"""


def infer_sequence_type(sequence: str) -> str:
    """依据字母表推断类型；不能确信时返回 nucleic 或 unknown，不硬猜 DNA/RNA。"""
    raw = PLACEHOLDER_RE.sub("", sequence)
    letters = {ch.upper() for ch in raw if ch.isalpha()}
    if not letters:
        return "unknown"
    if letters <= NUCLEIC_LETTERS:
        has_t, has_u = "T" in letters, "U" in letters
        if has_t and has_u:
            return "nucleic_mixed"
        if has_u:
            return "rna"
        if has_t:
            return "dna"
        return "nucleic"
    protein_letters = (PROTEIN_BASES | PROTEIN_AMBIG | PROTEIN_MODIFICATIONS) - {"*"}
    if letters <= protein_letters:
        return "protein"
    return "unknown"


def resolve_alphabet(seq_type: str, alphabet: str) -> str:
    """把 auto 解析为具体字母表档案。"""
    if alphabet != "auto":
        return alphabet
    return {
        "dna": "dna",
        "rna": "rna",
        "nucleic": "nucleic-mixed",
        "nucleic_mixed": "nucleic-mixed",
        "protein": "protein",
    }.get(seq_type, "nucleic-mixed")


def clean_sequence(sequence: str, seq_type: str = "auto",
                   placeholder_mode: str = "expand", alphabet: str = "auto",
                   gaps: str = "substitute", unknown: str = "substitute"):
    """把序列映射进目标字母表；返回 (结果, 映射审计, 类型, 字母表名)。

    保证：结果中每个字符都 ∈ 该字母表的 allowed（除非按 drop 删除）。
    """
    if placeholder_mode not in PLACEHOLDER_MODES:
        raise ValueError(f"未知占位符模式：{placeholder_mode}")
    if alphabet not in ALPHABET_CHOICES:
        raise ValueError(f"未知字母表：{alphabet}")
    if gaps not in GAP_MODES:
        raise ValueError(f"未知 gap 模式：{gaps}")
    if unknown not in UNKNOWN_MODES:
        raise ValueError(f"未知未知字符模式：{unknown}")

    actual_type = infer_sequence_type(sequence) if seq_type == "auto" else seq_type
    profile_name = resolve_alphabet(actual_type, alphabet)
    profile = ALPHABET_PROFILES[profile_name]
    allowed, fallback = profile["allowed"], profile["fallback"]

    stats = collections.Counter()
    out = []
    i = 0
    while i < len(sequence):
        match = PLACEHOLDER_RE.match(sequence, i)
        if match:
            k = int(match.group(1))
            if placeholder_mode == "expand":
                out.append(fallback * k)
                stats["占位符→fallback×k"] += 1
                stats["占位符还原字符数"] += k
            else:
                out.append(match.group(0))
                stats["占位符(保留·人工复核)"] += 1
            i = match.end()
            continue

        ch = sequence[i]
        upper = ch.upper()
        if upper in allowed:
            if ch.islower() and upper.isalpha():
                out.append(upper)
                stats["小写→大写"] += 1
            else:
                out.append(upper)
        elif ch in GAPS:
            if gaps == "drop":
                stats["gap(删除)"] += 1
            else:
                out.append(fallback)
                stats["gap→fallback"] += 1
        elif upper in IUPAC_AMBIG or upper in PROTEIN_AMBIG or upper in PROTEIN_MODIFICATIONS:
            # 字母表不含该模糊码/特殊状态
            if unknown == "drop":
                stats["未知(删除)"] += 1
            else:
                out.append(fallback)
                stats[f"{upper}→{fallback}"] += 1
        else:
            if unknown == "drop":
                stats["未知(删除)"] += 1
            else:
                out.append(fallback)
                stats[f"未知字符→{fallback}"] += 1
        i += 1

    cleaned = "".join(out)
    return cleaned, stats, actual_type, profile_name


def agent_clean(sequence: str, seq_type: str = "auto",
                placeholder_mode: str = "expand"):
    """向后兼容接口：返回 (结果, 审计)。"""
    cleaned, stats, _, _ = clean_sequence(sequence, seq_type, placeholder_mode)
    return cleaned, stats


def restore_aligned_column(original_seq: str, original_align: str,
                           placeholder_mode: str = "expand") -> str:
    """把与序列逐位对齐的注释列（如二级结构串）按同一占位符规则同步还原。

    注释列不是模型输入，因此只做长度对齐处理，不做字母表映射。
    """
    if placeholder_mode != "expand" or len(original_align) != len(original_seq):
        return original_align
    out = []
    i = 0
    while i < len(original_seq):
        m = PLACEHOLDER_RE.match(original_seq, i)
        if m:
            k = int(m.group(1))
            span = original_align[i:m.end()]
            out.append(span[0] * k if span and len(set(span)) == 1 else "." * k)
            i = m.end()
            continue
        out.append(original_align[i])
        i += 1
    return "".join(out)


def _sequence_profile(value: str):
    """列级检测：返回 (是否可信, 类型, 字母表覆盖率, 长度)。"""
    without_placeholder = PLACEHOLDER_RE.sub("", value)
    remaining = without_placeholder.replace("\r", "").replace("\n", "")
    if any(ch.isspace() for ch in remaining):
        return False, "unknown", 0.0, len(value)
    letters = [ch.upper() for ch in remaining if ch.isalpha()]
    if len(letters) < 10:
        return False, "unknown", 0.0, len(value)
    relevant = [ch for ch in remaining if ch not in GAPS]
    if not relevant:
        return False, "unknown", 0.0, len(value)
    nucleic_hits = sum((ch.upper() in NUCLEIC_LETTERS) if ch.isalpha() else False
                       for ch in relevant)
    protein_alphabet = PROTEIN_BASES | PROTEIN_AMBIG | PROTEIN_MODIFICATIONS
    protein_hits = sum((ch.upper() in protein_alphabet) if ch.isalpha() else ch == "*"
                       for ch in relevant)
    coverage = max(nucleic_hits, protein_hits) / len(relevant)
    seq_type = infer_sequence_type(value)
    return (seq_type != "unknown" and coverage >= 0.97), seq_type, coverage, len(value)


def discover_sequence_columns(fieldnames, sample_rows, minimum_average_length=20):
    """自动检测全部可信序列列；宁少勿滥，不能可靠识别时返回空列表。"""
    candidates = []
    minimum_values = min(5, max(1, len(sample_rows)))
    for col in fieldnames:
        values = []
        for row in sample_rows:
            value = row.get(col) or ""
            if value:
                values.append(value)
                if len(values) >= 200:
                    break
        if len(values) < minimum_values:
            continue
        profiles = [_sequence_profile(v) for v in values]
        good = [p for p in profiles if p[0]]
        if len(good) / len(values) < 0.8:
            continue
        average_length = sum(len(v) for v in values) / len(values)
        if average_length < minimum_average_length:
            continue
        type_counts = collections.Counter(p[1] for p in good)
        candidates.append({
            "column": col,
            "type_hint": type_counts.most_common(1)[0][0],
            "non_empty_samples": len(values),
            "plausible_samples": len(good),
            "average_length": round(average_length, 1),
            "alphabet_coverage": round(sum(p[2] for p in good) / len(good), 4),
            "detection": "auto_heuristic",
        })
    return candidates


def _parser():
    parser = argparse.ArgumentParser(
        description="把生物序列清洗成可直接进模型输入口的数据（保证输出字符全合法）")
    parser.add_argument("csv", help="输入 CSV 路径")
    parser.add_argument("--seq-col", dest="seq_cols", action="append", default=None,
                        help="只处理指定列；可重复使用，也支持逗号分隔")
    parser.add_argument("--seq-type", choices=TYPE_CHOICES, default="auto",
                        help="序列类型覆盖（默认逐条自动判断）")
    parser.add_argument("--alphabet", choices=ALPHABET_CHOICES, default="auto",
                        help="目标模型字母表（默认 auto：按类型选 dna/rna/protein）")
    parser.add_argument("--gaps", choices=GAP_MODES, default="substitute",
                        help="gap 处理：substitute 映射为 N/X（默认，保持对齐）；drop 删除")
    parser.add_argument("--unknown", choices=UNKNOWN_MODES, default="substitute",
                        help="其它非字母表字符：substitute 映射为 N/X（默认）；drop 删除")
    parser.add_argument("--placeholder", choices=PLACEHOLDER_MODES, default="expand",
                        help="占位符 *[k]*：expand 还原 k 个 N/X（默认，已实证）；keep 保留")
    parser.add_argument("--align-col", dest="align_cols", action="append", default=None,
                        help="需与序列逐位对齐的注释列（如 Structure）；可重复使用，"
                             "支持 序列列=注释列 显式配对")
    parser.add_argument("--out", default="curated_output", help="输出目录")
    parser.add_argument("--sample-rows", type=int, default=2000,
                        help="自动检测最多读取的开头行数，默认 2000")
    return parser


def run(args) -> dict:
    path = os.path.abspath(args.csv)
    if not os.path.isfile(path):
        raise CurateInputError(f"找不到文件：{path}")
    if args.sample_rows < 1:
        raise CurateInputError("--sample-rows 必须是正整数。")

    try:
        source = open(path, encoding="utf-8-sig", errors="replace", newline="")
    except OSError as exc:
        raise CurateInputError(f"无法读取输入文件：{exc}") from exc

    with source:
        reader = csv.DictReader(source)
        fieldnames = reader.fieldnames
        if not fieldnames:
            raise CurateInputError("CSV 没有有效表头，未写入任何结果。")
        sample_rows = []
        for i, row in enumerate(reader):
            if i >= args.sample_rows:
                break
            sample_rows.append(row)
    if not sample_rows:
        raise CurateInputError("CSV 只有表头、没有数据行，未写入任何结果。")

    if args.seq_cols:
        requested = [item.strip() for raw in args.seq_cols
                     for item in raw.split(",") if item.strip()]
        seq_cols = list(dict.fromkeys(requested))
        missing = [col for col in seq_cols if col not in fieldnames]
        if missing:
            available = "、".join(fieldnames[:20])
            more = f"……共 {len(fieldnames)} 列" if len(fieldnames) > 20 else ""
            raise CurateInputError(
                f"指定列不存在：{missing}。可用字段：{available}{more}。未写入任何结果。")
        detection = [{"column": col, "type_hint": args.seq_type,
                      "detection": "user_specified"} for col in seq_cols]
    else:
        detected = discover_sequence_columns(fieldnames, sample_rows)
        if not detected:
            available = "、".join(fieldnames[:20])
            more = f"……共 {len(fieldnames)} 列" if len(fieldnames) > 20 else ""
            raise CurateInputError(
                "未能可靠识别任何核酸/蛋白质序列列，已停止，未写入结果，避免误处理。"
                f"可用字段：{available}{more}。请用 --seq-col 列名 显式指定。")
        seq_cols = [item["column"] for item in detected]
        detection = detected

    align_cols, align_pairs = [], []
    if args.align_cols:
        requested_align = [item.strip() for raw in args.align_cols
                           for item in raw.split(",") if item.strip()]
        for spec in dict.fromkeys(requested_align):
            if "=" in spec:
                seq_part, align_part = (p.strip() for p in spec.split("=", 1))
                if seq_part not in seq_cols or align_part not in fieldnames:
                    raise CurateInputError(
                        f"--align-col 配对无效：{spec}。左侧须是已处理序列列 {seq_cols}，"
                        "右侧须是存在的列。")
                align_cols.append(align_part)
                align_pairs.append((seq_part, align_part))
                continue
            if spec not in fieldnames:
                raise CurateInputError(f"--align-col 指定列不存在：{spec}。未写入任何结果。")
            best, best_rate = None, 0.0
            for scol in seq_cols:
                both = equal = 0
                for row in sample_rows:
                    a, b = row.get(scol) or "", row.get(spec) or ""
                    if a and b:
                        both += 1
                        equal += int(len(a) == len(b))
                rate = equal / both if both else 0.0
                if rate > best_rate:
                    best, best_rate = scol, rate
            if best is None or best_rate < 0.5:
                raise CurateInputError(
                    f"--align-col {spec} 无法与任何序列列自动配对"
                    f"（最高原长一致率 {best_rate:.0%}）。请用 --align-col 序列列={spec}。")
            align_cols.append(spec)
            align_pairs.append((best, spec))
        align_cols = list(dict.fromkeys(align_cols))
    align_for_seq = collections.defaultdict(list)
    for scol, acol in align_pairs:
        align_for_seq[scol].append(acol)

    os.makedirs(args.out, exist_ok=True)
    out_csv = os.path.join(args.out, "cleaned.csv")
    out_summary = os.path.join(args.out, "curation_summary.json")
    out_fields = list(fieldnames)
    for col in seq_cols:
        out_fields.extend([f"{col}_curated", f"{col}_seq_type",
                           f"{col}_alphabet", f"{col}_curation_summary"])
    for col in align_cols:
        out_fields.append(f"{col}_curated")
    out_fields.append("curation_summary")

    overall = collections.Counter()
    per_column = {col: collections.Counter() for col in seq_cols}
    per_column_types = {col: collections.Counter() for col in seq_cols}
    per_column_alpha = {col: collections.Counter() for col in seq_cols}
    per_column_nonempty = collections.Counter()
    allowed_union = {col: set() for col in seq_cols}
    n_rows = 0
    lengths_changed = 0
    placeholders_seen = 0
    n_restored = 0
    align_misaligned = 0
    in_loop_violations = 0

    with open(path, encoding="utf-8-sig", errors="replace", newline="") as src, \
            open(out_csv, "w", encoding="utf-8-sig", newline="") as dst:
        reader = csv.DictReader(src)
        writer = csv.DictWriter(dst, fieldnames=out_fields, extrasaction="ignore")
        writer.writeheader()
        for row in reader:
            n_rows += 1
            row_details = []
            for col in seq_cols:
                value = row.get(col) or ""
                if not value:
                    for suffix in ("_curated", "_seq_type", "_alphabet",
                                   "_curation_summary"):
                        row[f"{col}{suffix}"] = ""
                    continue
                cleaned, stats, actual_type, profile_name = clean_sequence(
                    value, args.seq_type, args.placeholder, args.alphabet,
                    args.gaps, args.unknown)
                allowed = ALPHABET_PROFILES[profile_name]["allowed"]
                allowed_union[col] |= allowed
                # 即时保证：输出字符必须全在字母表内
                bad = {c for c in cleaned if c not in allowed}
                if bad:
                    in_loop_violations += 1
                n_ph = (stats.get("占位符→fallback×k", 0)
                        + stats.get("占位符(保留·人工复核)", 0))
                placeholders_seen += n_ph
                n_restored += stats.get("占位符还原字符数", 0)
                if len(cleaned) != len(value):
                    lengths_changed += 1
                overall.update(stats)
                per_column[col].update(stats)
                per_column_types[col][actual_type] += 1
                per_column_alpha[col][profile_name] += 1
                per_column_nonempty[col] += 1
                detail = {
                    "column": col,
                    "sequence_type": actual_type,
                    "target_alphabet": profile_name,
                    "original_length": len(value),
                    "curated_length": len(cleaned),
                    "mapping": dict(stats),
                }
                if n_ph:
                    detail["placeholders"] = n_ph
                row[f"{col}_curated"] = cleaned
                row[f"{col}_seq_type"] = actual_type
                row[f"{col}_alphabet"] = profile_name
                row_details.append(detail)

                for acol in align_for_seq.get(col, []):
                    avalue = row.get(acol) or ""
                    if not avalue:
                        row[f"{acol}_curated"] = ""
                        continue
                    if len(avalue) != len(value):
                        align_misaligned += 1
                        row[f"{acol}_curated"] = avalue
                        detail.setdefault("align_columns", {})[acol] = "跳过：原长不一致"
                        continue
                    fixed = restore_aligned_column(value, avalue, args.placeholder)
                    row[f"{acol}_curated"] = fixed
                    detail.setdefault("align_columns", {})[acol] = {
                        "original_length": len(avalue), "curated_length": len(fixed)}
                row[f"{col}_curation_summary"] = json.dumps(detail, ensure_ascii=False)
            row["curation_summary"] = json.dumps(
                {"columns": row_details}, ensure_ascii=False) if row_details else ""
            writer.writerow(row)

    # 独立复核：重新读回成品，逐字符确认全部落在目标字母表内
    violations = collections.Counter()
    checked_chars = 0
    with open(out_csv, encoding="utf-8-sig", newline="") as check:
        for out_row in csv.DictReader(check):
            for col in seq_cols:
                value = out_row.get(f"{col}_curated") or ""
                if not value:
                    continue
                allowed = allowed_union[col]
                checked_chars += len(value)
                for ch in value:
                    if ch not in allowed:
                        violations[f"{col}:{ch}"] += 1

    summary = {
        "input": path,
        "rows": n_rows,
        "sequence_columns": seq_cols,
        "seq_type_override": args.seq_type,
        "target_alphabet": args.alphabet,
        "gap_mode": args.gaps,
        "unknown_mode": args.unknown,
        "placeholder_mode": args.placeholder,
        "align_columns": align_cols,
        "align_pairs": [{"sequence_column": s, "align_column": a} for s, a in align_pairs],
        "column_detection": detection,
        "mapping": dict(overall),
        "mapping_reasons": DECISION_REASONS,
        "columns": {
            col: {
                "non_empty_rows": per_column_nonempty[col],
                "types": dict(per_column_types[col]),
                "alphabets": dict(per_column_alpha[col]),
                "allowed_chars": "".join(sorted(allowed_union[col])),
                "mapping": dict(per_column[col]),
            } for col in seq_cols
        },
        "placeholders_seen": placeholders_seen,
        "chars_restored": n_restored,
        "rows_with_length_change": lengths_changed,
        "align_column_skipped": align_misaligned,
        "alphabet_check": {
            "checked_chars": checked_chars,
            "violations": dict(violations),
            "passed": not violations and in_loop_violations == 0,
        },
        "notes": [
            "输出保证只含目标字母表内的字符：非标准字符被映射为 N/X 或按选项删除。",
            "IUPAC 模糊码在严格字母表下映射为 N/X（与 cleanfa 选择一致）；"
            "若 tokenizer 支持，用 --alphabet dna-iupac / rna-iupac 原样保留。",
            "gap 默认映射为 N/X 以保持与结构注释的逐位对齐；--gaps drop 则删除。",
            "占位符 *[k]* 按 k 个 N/X 还原（k 已用未损坏 FASTA 侧翼锚定实证）。",
            "小写代表低置信度/可变位点，转大写后是该字母的合法 token；映射计入审计。",
            "每行每列的映射明细写在 <列名>_curation_summary，可逐字符复查。",
            "注释列（--align-col，如 Structure）不是模型输入，只做长度对齐，不做字母表映射。",
            "原文件不改动。",
        ],
    }
    with open(out_summary, "w", encoding="utf-8") as dst:
        json.dump(summary, dst, ensure_ascii=False, indent=2)
    summary["out_csv"] = out_csv
    summary["out_summary"] = out_summary
    return summary


def main(argv=None):
    args = _parser().parse_args(argv)
    try:
        summary = run(args)
    except CurateInputError as exc:
        print(f"[错误] {exc}")
        return 2
    except (OSError, csv.Error) as exc:
        print(f"[错误] 读写 CSV 失败：{exc}；请检查输入格式与输出目录权限。")
        return 2

    print("=" * 76)
    print("清理完成：输出可直接进模型输入口（原文件未改动）")
    print("=" * 76)
    print(f"输入    : {summary['input']}")
    print(f"行数    : {summary['rows']:,}")
    print(f"序列列  : {len(summary['sequence_columns'])} 个"
          f"（目标字母表 {summary['target_alphabet']}，gap {summary['gap_mode']}）")
    for col in summary["sequence_columns"]:
        d = summary["columns"][col]
        print(f"  - {col}: {d['non_empty_rows']:,} 条，类型 {d['types']}，"
              f"字母表 {d['alphabets']}")
        print(f"      允许字符集: {d['allowed_chars']}")
    print("映射汇总（原来 → 变成什么）:")
    for key, value in sorted(summary["mapping"].items(), key=lambda x: -x[1]):
        print(f"  {value:>9,}  {key}")
    if summary["placeholders_seen"]:
        print(f"占位符 {summary['placeholders_seen']:,} 个 → 还原 "
              f"{summary['chars_restored']:,} 个字符；"
              f"长度变化行 {summary['rows_with_length_change']:,}")
    check = summary["alphabet_check"]
    print(f"字母表复核: 检查 {check['checked_chars']:,} 个字符，"
          f"非法 {sum(check['violations'].values()):,}"
          f" → {'通过 ✅' if check['passed'] else '未通过 ❌ ' + str(check['violations'])}")
    print(f"成品    : {summary['out_csv']}")
    print(f"摘要    : {summary['out_summary']}")
    return 0 if check["passed"] else 3


if __name__ == "__main__":
    sys.exit(main())
