"""未知数据约定的发现模块 —— 这是本仓库唯一真正需要 LLM 的地方。

为什么需要 LLM：
    清洗本身（把字符映射进词表）不需要模型，规则更快更稳更便宜——这一点本仓库从不含糊。
    真正需要模型的是**面对一份陌生的脏数据，判断它的记号是什么意思**。
    脚本只能执行人已经想好的规则；没人告诉它 `*[k]*` 是什么时，它无能为力。

但也不能让模型说了算（会幻觉）。所以这里定下一个不可让步的分工：

    LLM 负责"提出可证伪的假设"  →  代码负责"跑实验判定真假"

    ┌─ 阶段 1  观察（确定性）：扫描数据，找出可疑记号及其分布
    ├─ 阶段 2  提议（LLM）  ：对每个记号提出假设，并从**受支持预测清单**中选一个
    ├─ 阶段 3  验证（确定性）：用参照列做锚定等实验，量出真实结果
    └─ 阶段 4  裁定（确定性）：预测成立 → CONFIRMED；不成立 → REFUTED；
                              无法验证 → UNTESTABLE（绝不当成通过）

LLM 说不出"用什么实验验证"时，该假设直接标为 UNTESTABLE —— 不允许"看起来像"就下结论。

用法：
    # 离线自检（内置 mock 假设，不连任何模型）
    python discover_conventions.py <raw.csv> --seq-col Sequence --ref-col FASTA_sequence --llm mock

    # 接入任意模型：命令从 stdin 读 prompt、把回复写到 stdout
    python discover_conventions.py <raw.csv> --seq-col Sequence --ref-col FASTA_sequence \
        --llm-cmd "python my_llm_bridge.py"

    # 把确认的结论输出成规则建议，人工确认后再给 curate.py 用
    python discover_conventions.py <raw.csv> ... --emit-rules rules.json
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import os
import re
import subprocess
import sys

# ---------------------------------------------------------------------------
# 支持的可证伪预测：LLM 只能从这个清单里选，避免它编造无法执行的验证方案
# ---------------------------------------------------------------------------
SUPPORTED_PREDICTIONS = {
    "span_equals_token_width":
        "锚定出的真实跨度 == 该记号的文本宽度（说明记号只是注解，不占真实碱基）",
    "span_equals_number_in_token":
        "锚定出的真实跨度 == 记号里的数字（说明数字是被替换掉的碱基数）",
    "span_equals_one":
        "锚定出的真实跨度 == 1（说明该处恰好对应一个真实碱基）",
    "span_is_zero":
        "锚定出的真实跨度 == 0（说明该处没有真实碱基）",
    "span_exceeds_token_width":
        "真实跨度 > 记号文本宽度（说明记号展开后应比自身更长）",
}

CONFIRM_THRESHOLD = 0.95
REFUTE_THRESHOLD = 0.05
MIN_SAMPLES = 20          # 少于这个数不足以判定，标 UNTESTABLE


# ---------------------------------------------------------------------------
# 阶段 1：观察（确定性，不需要模型）
# ---------------------------------------------------------------------------
SUSPECT_PATTERNS = [
    ("star_bracket_number", re.compile(r"\*\[ *(\d+)\]\*"),
     "星号方括号数字记号，如 *[19]*"),
    ("lowercase_letters", re.compile(r"[a-z]"),
     "小写字母（比对惯例中常代表低置信度/可变位点）"),
    ("gap_characters", re.compile(r"[-.]"),
     "比对缺口字符 - 或 ."),
    ("digit_run", re.compile(r"\d+"),
     "裸露的数字串"),
    ("ambiguous_code", re.compile(r"[NRYWSKMBDHVX]"),
     "IUPAC 模糊码"),
    ("non_alnum", re.compile(r"[^A-Za-z0-9]"),
     "其它非字母数字符号"),
]


def observe(csv_path, seq_col, sample_limit):
    """扫描数据，返回 {列: {字符直方图, 命中的可疑记号}}。"""
    profile = {}
    with open(csv_path, encoding="utf-8-sig", errors="replace", newline="") as f:
        reader = csv.DictReader(f)
        if seq_col not in (reader.fieldnames or []):
            raise SystemExit(f"[错误] 列不存在：{seq_col}")
        hist = collections.Counter()
        found = collections.OrderedDict()
        examples = {}
        rows = 0
        for row in reader:
            value = row.get(seq_col) or ""
            if not value:
                continue
            rows += 1
            hist.update(value)
            for name, regex, desc in SUSPECT_PATTERNS:
                matches = regex.findall(value) if regex.groups else regex.findall(value)
                n = len(matches)
                if not n:
                    continue
                if name not in found:
                    found[name] = {"pattern": name, "description": desc,
                                   "regex": regex.pattern, "count": 0, "distinct": set()}
                found[name]["count"] += n
                for m in matches[:5]:
                    token = m if isinstance(m, str) else m[0]
                    found[name]["distinct"].add(token[:12])
                if name not in examples:
                    hit = regex.search(value)
                    lo = max(0, hit.start() - 25)
                    examples[name] = value[lo:hit.end() + 25]
            if rows >= sample_limit:
                break

    suspects = []
    for name, info in found.items():
        info["distinct"] = sorted(info["distinct"])[:8]
        info["example"] = examples.get(name, "")
        suspects.append(info)
    profile[seq_col] = {
        "rows_scanned": rows,
        "alphabet": "".join(sorted(hist)),
        "char_histogram": dict(hist.most_common(15)),
        "suspects": suspects,
    }
    return profile


# ---------------------------------------------------------------------------
# 阶段 2：提议（LLM）
# ---------------------------------------------------------------------------
PROMPT_TEMPLATE = """你是资深生物信息数据管家。下面是一份序列数据里自动扫描出的可疑记号。
请对每个记号提出一个关于"它到底是什么"的假设，并指定用哪个预测来验证它。

只能使用下列预测之一（必须原样照抄字串）：
{predictions}

可用的验证手段：
- 若提供了参照列（未被污染的同类序列），将用侧翼锚定量出记号的真实跨度，再核对上面选定的预测。

要求：
1. 每个记号给一条假设；不确定就给多条不同假设（会被分别验证）。
2. 不许编造预测名；不许给无法验证的假设（那种情况请把 prediction 留空字符串 ""）。
3. rationale 说明你为什么这么猜，一句话。

只输出 JSON 数组，元素格式：
{{"pattern": "<记号名>", "hypothesis": "<计数|固定宽度注解|低置信度位点|修饰|噪声|其它>",
  "rationale": "<一句话>", "prediction": "<上面清单里的字串，或空串>"}}

扫描结果：
{profile}
"""


def _mock_llm(_prompt: str) -> str:
    """离线自检用：返回一组固定的、覆盖正反面的假设。"""
    return json.dumps([
        {"pattern": "star_bracket_number", "hypothesis": "计数",
         "rationale": "数字很像被替换掉的碱基个数", "prediction": "span_equals_number_in_token"},
        {"pattern": "star_bracket_number", "hypothesis": "固定宽度注解",
         "rationale": "它也可能只是个不占真实位置的标记", "prediction": "span_equals_token_width"},
        {"pattern": "lowercase_letters", "hypothesis": "低置信度位点",
         "rationale": "比对惯例中小写是低置信度", "prediction": "span_equals_token_width"},
    ], ensure_ascii=False)


def call_llm(mode, cmd, prompt):
    if mode == "mock":
        return _mock_llm(prompt)
    if not cmd:
        raise SystemExit("[错误] 需要 --llm mock 或 --llm-cmd '<命令>'")
    try:
        result = subprocess.run(cmd, shell=True, input=prompt, capture_output=True,
                                text=True, encoding="utf-8", errors="replace")
    except OSError as exc:
        raise SystemExit(f"[错误] 调用模型失败：{exc}")
    if result.returncode != 0:
        raise SystemExit(f"[错误] 模型命令返回 {result.returncode}：{result.stderr[:400]}")
    return result.stdout


def parse_hypotheses(raw):
    text = raw.strip()
    if text.startswith("```"):
        text = re.sub(r"^```[a-zA-Z]*\n?|```$", "", text).strip()
    try:
        data = json.loads(text)
    except json.JSONDecodeError:
        return [], "模型输出不是合法 JSON"
    if not isinstance(data, list):
        return [], "模型输出不是数组"
    cleaned = []
    for item in data:
        if not isinstance(item, dict):
            continue
        pred = item.get("prediction") or ""
        if pred and pred not in SUPPORTED_PREDICTIONS:
            item["prediction"] = ""
            item["rationale"] = (item.get("rationale", "") +
                                 " [原预测名不受支持，已置空]").strip()
        cleaned.append(item)
    return cleaned, ""


# ---------------------------------------------------------------------------
# 阶段 3：验证（确定性探测器）
# ---------------------------------------------------------------------------
def rna_to_dna(text):
    return text.replace("U", "T").replace("u", "t")


LETTER_RE = re.compile(r"[A-Za-z]")
MIN_FLANK = 8   # 侧翼至少要有这么多字母才够定位


def _flank(seq, start, end, direction, want):
    """取目标记号旁的侧翼，**只收字母**并在收够 want 个后停止。

    只收字母是刻意的：参照列是干净的序列，只有字母能在里面定位；
    一旦把 `*`、`[`、`-` 这些噪声收进侧翼，锚点就永远找不到——
    这会静默地把样本数变成 0（本模块早期的真实缺陷）。
    """
    out = []
    i = start - 1 if direction < 0 else end
    while 0 <= i < len(seq) and len(out) < want:
        if LETTER_RE.match(seq[i]):
            out.append(seq[i])
            i += direction
        else:
            break          # 碰到噪声就停，宁可侧翼短也不污染锚点
    if len(out) < MIN_FLANK:
        return None
    if direction < 0:
        out.reverse()      # 向左收集是逆序的，必须翻回来才能当锚点
    return rna_to_dna("".join(out)).upper()


def probe_anchor_span(csv_path, seq_col, ref_col, regex_src, flank, sample_limit):
    """用参照列做侧翼锚定，量出每个记号的真实跨度。

    返回 {"tested": n, "spans": [...], "records": [(token, span), ...]}
    records 同时用于需要 token 的预测（宽度、数字），因此一次扫描即可。
    """
    regex = re.compile(regex_src)
    spans = []
    records = []
    with open(csv_path, encoding="utf-8-sig", errors="replace", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            seq = row.get(seq_col) or ""
            ref = (row.get(ref_col) or "").upper()
            if not seq or not ref:
                continue
            for m in regex.finditer(seq):
                left = _flank(seq, m.start(), m.end(), -1, flank)
                right = _flank(seq, m.start(), m.end(), +1, flank)
                if left is None or right is None:
                    continue
                li = ref.find(left)
                if li < 0:
                    continue
                ri = ref.find(right, li + len(left))
                if ri < 0:
                    continue
                span = ri - (li + len(left))
                spans.append(span)
                records.append((m.group(0), span))
            if len(spans) >= sample_limit:
                break
    return {"tested": len(spans), "spans": spans, "records": records}


def evaluate_predictions(probe_result):
    """对每个受支持的预测算出命中率。返回 {预测名: 命中率 或 None}。"""
    spans = probe_result["spans"]
    records = probe_result["records"]
    if not spans:
        return {}
    n = len(spans)
    out = {}

    hit_w = sum(1 for token, span in records if span == len(token))
    out["span_equals_token_width"] = hit_w / len(records)

    numeric = [(token, span) for token, span in records if re.findall(r"\d+", token)]
    if numeric:
        hit_num = sum(1 for token, span in numeric
                      if span == int(re.findall(r"\d+", token)[0]))
        out["span_equals_number_in_token"] = hit_num / len(numeric)
    else:
        out["span_equals_number_in_token"] = None

    out["span_equals_one"] = sum(1 for s in spans if s == 1) / n
    out["span_is_zero"] = sum(1 for s in spans if s == 0) / n
    out["span_exceeds_token_width"] = (
        sum(1 for token, span in records if span > len(token)) / len(records))
    return out


def verdict(rate, tested):
    if tested < MIN_SAMPLES:
        return "UNTESTABLE"
    if rate is None:
        return "UNTESTABLE"
    if rate >= CONFIRM_THRESHOLD:
        return "CONFIRMED"
    if rate <= REFUTE_THRESHOLD:
        return "REFUTED"
    return "WEAK"


# ---------------------------------------------------------------------------
# 阶段 4：裁定与报告
# ---------------------------------------------------------------------------
def main(argv=None):
    parser = argparse.ArgumentParser(
        description="发现未知数据约定的含义（LLM 提假设，代码做验证）")
    parser.add_argument("csv")
    parser.add_argument("--seq-col", required=True)
    parser.add_argument("--ref-col", default=None,
                        help="未被污染的参照列；没有它则锚定实验无法进行")
    parser.add_argument("--llm", choices=["mock"], default=None,
                        help="用内置 mock 假设离线自检")
    parser.add_argument("--llm-cmd", default=None,
                        help="模型桥命令：从 stdin 读 prompt、回复写 stdout")
    parser.add_argument("--flank", type=int, default=15)
    parser.add_argument("--sample-rows", type=int, default=20000)
    parser.add_argument("--anchor-limit", type=int, default=3000,
                        help="锚定实验最多测多少个记号")
    parser.add_argument("--emit-rules", default=None,
                        help="把确认的结论写成规则建议 JSON")
    args = parser.parse_args(argv)

    if not args.llm and not args.llm_cmd:
        print("[错误] 请指定 --llm mock 或 --llm-cmd '<命令>'")
        return 2

    print("阶段 1/4  观察数据（确定性扫描）……")
    profile = observe(args.csv, args.seq_col, args.sample_rows)
    col_profile = profile[args.seq_col]
    print(f"  扫描 {col_profile['rows_scanned']:,} 行")
    print(f"  字符集: {col_profile['alphabet'][:60]}")
    for s in col_profile["suspects"]:
        print(f"  - {s['pattern']}: {s['count']:,} 次  例: {s['example'][:48]}")

    print("\n阶段 2/4  请模型提出可证伪的假设……")
    prompt = PROMPT_TEMPLATE.format(
        predictions="\n".join(f"  - {k}：{v}" for k, v in SUPPORTED_PREDICTIONS.items()),
        profile=json.dumps(col_profile, ensure_ascii=False, indent=1),
    )
    raw = call_llm(args.llm, args.llm_cmd, prompt)
    hypotheses, err = parse_hypotheses(raw)
    if err:
        print(f"  [警告] {err} → 无法提出假设，结束（不猜）")
        return 2
    print(f"  收到 {len(hypotheses)} 条假设")

    print("\n阶段 3/4  跑验证实验（确定性，模型不参与判定）……")
    by_pattern = {s["pattern"]: s for s in col_profile["suspects"]}
    findings = []
    for h in hypotheses:
        pattern = h.get("pattern", "")
        info = by_pattern.get(pattern)
        prediction = h.get("prediction", "")
        entry = {"pattern": pattern, "hypothesis": h.get("hypothesis", ""),
                 "rationale": h.get("rationale", ""), "prediction": prediction}
        if info is None:
            entry.update(verdict="UNTESTABLE", reason="该记号在数据中不存在", tested=0)
            findings.append(entry); continue
        if not prediction:
            entry.update(verdict="UNTESTABLE", reason="模型未给出可验证的预测", tested=0)
            findings.append(entry); continue
        if not args.ref_col:
            entry.update(verdict="UNTESTABLE", reason="未提供参照列，无法锚定", tested=0)
            findings.append(entry); continue

        probe = probe_anchor_span(args.csv, args.seq_col, args.ref_col,
                                  info["regex"], args.flank, args.anchor_limit)
        rates = evaluate_predictions(probe)
        rate = rates.get(prediction)
        v = verdict(rate, probe["tested"])
        entry.update(
            verdict=v, tested=probe["tested"],
            match_rate=(round(rate, 4) if rate is not None else None),
            span_examples=probe["records"][:3],
            all_rates={k: (round(x, 4) if x is not None else None)
                       for k, x in rates.items()},
        )
        findings.append(entry)

    print("\n阶段 4/4  裁定")
    print("=" * 84)
    print(f"{'记号':<22}{'假设':<14}{'预测':<28}{'样本':>7}{'命中率':>9}{'裁定':>12}")
    print("-" * 84)
    for f in findings:
        rate = f.get("match_rate")
        print(f"{f['pattern']:<22}{f.get('hypothesis','')[:12]:<14}"
              f"{f.get('prediction','')[:26]:<28}{f.get('tested',0):>7}"
              f"{(f'{100*rate:.1f}%' if rate is not None else '—'):>9}"
              f"{f['verdict']:>12}")
    print("-" * 84)
    print("CONFIRMED = 预测成立；REFUTED = 预测被证伪；"
          f"WEAK = 证据不足；UNTESTABLE = 无法验证（不足 {MIN_SAMPLES} 样本或没有实验方案）")
    print()
    for f in findings:
        if f.get("reason"):
            print(f"  · {f['pattern']}/{f.get('hypothesis','')}: {f['verdict']} —— {f['reason']}")
    for f in findings:
        if f["verdict"] in ("CONFIRMED", "REFUTED") and f.get("all_rates"):
            print(f"  · {f['pattern']} 各预测命中率: {f['all_rates']}")

    if args.emit_rules:
        rules = _suggest_rules(findings)
        with open(args.emit_rules, "w", encoding="utf-8") as fh:
            json.dump({"source": args.csv, "sequence_column": args.seq_col,
                       "reference_column": args.ref_col,
                       "findings": findings, "suggested_rules": rules,
                       "note": "建议需人工确认后再交给 curate.py；未经验证的假设不得生效。"},
                      fh, ensure_ascii=False, indent=2)
        print(f"\n规则建议已写入: {args.emit_rules}（需人工确认）")
    return 0


def _suggest_rules(findings):
    """把 CONFIRMED 的结论翻译成 curate.py 可选用的规则建议。"""
    rules = []
    for f in findings:
        if f["verdict"] != "CONFIRMED":
            continue
        pred = f.get("prediction")
        if f["pattern"] == "star_bracket_number" and pred == "span_equals_number_in_token":
            rules.append({"regex": r"\*\[ *(\d+)\]\*", "action": "expand_to_n_by_number",
                          "evidence": f"锚定 {f['tested']} 例，命中率 {f['match_rate']}",
                          "maps_to": "占位符 -> k 个 N"})
        elif pred == "span_equals_token_width" and f["pattern"] == "lowercase_letters":
            rules.append({"regex": "[a-z]", "action": "uppercase_keep",
                          "evidence": f"锚定 {f['tested']} 例，命中率 {f['match_rate']}",
                          "maps_to": "小写 -> 同一碱基的大写（非插入记号）"})
    return rules


if __name__ == "__main__":
    sys.exit(main())
