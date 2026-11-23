"""curate.py 入口回归测试。

重点验证**目标**：输出必须全部落在目标字母表内（模型喂得进去），
同时映射可审计、失败不静默。使用合成小 CSV，不依赖真实用户数据。
"""
import contextlib
import csv
import io
import json
import os
import tempfile

from curate import (
    ALPHABET_PROFILES,
    CurateInputError,
    _parser,
    clean_sequence,
    discover_sequence_columns,
    infer_sequence_type,
    main,
    restore_aligned_column,
    run,
)

FAILURES = []


def check(name, condition):
    print(("PASS  " if condition else "FAIL  ") + name)
    if not condition:
        FAILURES.append(name)


def all_in(seq, alphabet):
    return set(seq) <= ALPHABET_PROFILES[alphabet]["allowed"]


# ---------------------------------------------------------------- 类型推断
check("DNA → dna", infer_sequence_type("ACGTNacgt") == "dna")
check("RNA → rna", infer_sequence_type("ACGUNacgu") == "rna")
check("无 T/U 不强行区分 DNA/RNA", infer_sequence_type("ACGRN") == "nucleic")
check("T/U 混合不互相视未知", infer_sequence_type("ACGTU") == "nucleic_mixed")
check("蛋白 U/O 特殊状态", infer_sequence_type("MKTFUISVAKBOXZUO") == "protein")

# ------------------------------------------------- 核心：输出必须合法
mixed_clean, mixed_stats, mixed_type, mixed_alpha = clean_sequence("ACGTUacgtu")
check("混合核酸小写→大写且全部合法",
      mixed_clean == "ACGTUACGTU" and all_in(mixed_clean, mixed_alpha))

prot_raw = "MKTNuobxj*.-"          # 含 U O B X J * gap
prot_clean, prot_stats, prot_type, prot_alpha = clean_sequence(prot_raw)
check("蛋白：非词表字符全部映射为 X，输出 100% 合法",
      prot_type == "protein" and all_in(prot_clean, prot_alpha)
      and prot_clean == "MKTNXXXXXXXX")
check("蛋白映射逐字符可审计",
      prot_stats["U→X"] == 1 and prot_stats["O→X"] == 1 and prot_stats["B→X"] == 1
      and prot_stats["J→X"] == 1 and prot_stats["*→X"] == 1
      and prot_stats["gap→fallback"] == 2)

dna_raw = "ACGTNacgt-*[ 2]*"       # 含小写 gap 占位符
dna_clean, dna_stats, dna_type, dna_alpha = clean_sequence(dna_raw)
check("DNA：gap 与占位符都变成合法 token",
      dna_type == "dna" and all_in(dna_clean, dna_alpha)
      and dna_clean == "ACGTNACGTNNN")
check("占位符按计数还原（k=2 → NN）而非保留标记",
      dna_stats["占位符→fallback×k"] == 1 and "*[" not in dna_clean)
check("gap 被映射而非保留", dna_stats["gap→fallback"] == 1 and "-" not in dna_clean)

# 严格字母表下模糊码必须映射为 N
amb_clean, amb_stats, _, amb_alpha = clean_sequence("ACGTRYWSNACGT")
check("严格字母表下 IUPAC 模糊码 → N",
      amb_clean == "ACGTNNNNNACGT" and all_in(amb_clean, amb_alpha)
      and amb_stats["R→N"] == 1 and amb_stats["Y→N"] == 1 and amb_stats["W→N"] == 1)

# tokenizer 支持模糊码时可原样保留
keep_clean, keep_stats, _, keep_alpha = clean_sequence(
    "ACGTRYWSNACGT", alphabet="dna-iupac")
check("--alphabet dna-iupac 时模糊码合法保留（不被映射）",
      keep_clean == "ACGTRYWSNACGT" and all_in(keep_clean, keep_alpha))

# 删除模式
drop_clean, drop_stats, _, drop_alpha = clean_sequence(
    "ACGT-ACGT", gaps="drop")
check("--gaps drop 时 gap 被删除且输出合法",
      drop_clean == "ACGTACGT" and all_in(drop_clean, drop_alpha)
      and drop_stats["gap(删除)"] == 1)
drop2_clean, _, _, drop2_alpha = clean_sequence("ACGT!ACGT", unknown="drop")
check("--unknown drop 时未知字符被删除且输出合法",
      drop2_clean == "ACGTACGT" and all_in(drop2_clean, drop2_alpha))

# 边界：k=0 / keep 模式 / 未知字母表
zero_clean, zero_stats, _, _ = clean_sequence("ACGU*[ 0]*ACGU")
check("k=0 时占位符消失", zero_clean == "ACGUACGU"
      and zero_stats["占位符还原字符数"] == 0)
kept_clean, kept_stats, _, _ = clean_sequence(dna_raw, placeholder_mode="keep")
check("keep 模式保留标记原样（供人工核对）",
      "*[ 2]*" in kept_clean and kept_stats["占位符(保留·人工复核)"] == 1)
try:
    clean_sequence("ACGT", alphabet="not-an-alphabet")
    bad_alpha = False
except ValueError:
    bad_alpha = True
check("未知字母表报错", bad_alpha)

# 对齐注释列
check("对齐注释列按同一计数还原",
      restore_aligned_column("ACGU*[19]*ACGU", "ACGU~~~~~~ACGU")
      == "ACGU" + "~" * 19 + "ACGU")
check("原长不一致时对齐列不冒险映射",
      restore_aligned_column("ACGU*[19]*ACGU", "SHORT") == "SHORT")

# ---------------------------------------------------------- 列自动检测
_RNA_SEQ = "ACGUNACGUNACGUNACGUNacgu*[ 2]*"
_STRUCT_PREFIX = "<<((....))>>"
sample = [{
    "id": f"item-{i}",
    "dna_seq": "ACGTNACGTNACGTNACGTNacgt-",
    "rna_seq": _RNA_SEQ,
    "protein_seq": "MKTFUISVAKBOXZUOMKNP*",
    "note": "this is not sequence data for machine learning",
    # 结构列必须与 rna_seq 等长，才能被自动配对并同步还原
    "structure": _STRUCT_PREFIX + "." * (len(_RNA_SEQ) - len(_STRUCT_PREFIX)),
} for i in range(6)]
fields = list(sample[0])
check("结构列与 rna_seq 原始等长", len(sample[0]["structure"]) == len(sample[0]["rna_seq"]))
candidates = discover_sequence_columns(fields, sample)
candidate_names = {item["column"] for item in candidates}
check("自动检测全部三种序列列",
      candidate_names == {"dna_seq", "rna_seq", "protein_seq"})
check("不误识别自然语言和结构字符串",
      "note" not in candidate_names and "structure" not in candidate_names
      and "id" not in candidate_names)

# ------------------------------------------------------------- 端到端
with tempfile.TemporaryDirectory() as tmp:
    source = os.path.join(tmp, "mini.csv")
    with open(source, "w", encoding="utf-8-sig", newline="") as dst:
        writer = csv.DictWriter(dst, fieldnames=fields)
        writer.writeheader()
        writer.writerows(sample)
    out_dir = os.path.join(tmp, "out")
    stdout = io.StringIO()
    with contextlib.redirect_stdout(stdout):
        code = main([source, "--out", out_dir, "--align-col", "structure"])
    check("入口成功返回 0", code == 0)
    check("报告声明通过字母表复核", "通过 ✅" in stdout.getvalue())
    with open(os.path.join(out_dir, "curation_summary.json"), encoding="utf-8") as src:
        summary = json.load(src)
    check("摘要记录全部三列", set(summary["sequence_columns"]) == candidate_names)
    check("字母表复核通过且检查了字符",
          summary["alphabet_check"]["passed"]
          and summary["alphabet_check"]["checked_chars"] > 0)
    check("占位符已还原并记录",
          summary["placeholders_seen"] > 0 and summary["chars_restored"] > 0)

    with open(os.path.join(out_dir, "cleaned.csv"), encoding="utf-8-sig", newline="") as src:
        rows = list(csv.DictReader(src))
    row = rows[0]
    # 全表扫描：任何序列列的产物都不得含字母表外字符
    leaks = []
    for r in rows:
        for col in candidate_names:
            value = r.get(f"{col}_curated") or ""
            allowed = ALPHABET_PROFILES[r.get(f"{col}_alphabet")]["allowed"]
            leaks += [c for c in value if c not in allowed]
    check("全表输出无任何字母表外字符（核心目标）", not leaks)
    check("gap 已消失", "-" not in (row["dna_seq_curated"] or ""))
    check("模糊码/特殊状态已消失",
          not (set(row["protein_seq_curated"]) & set("BJZUO*")))
    check("RNA 占位符已还原为 NN",
          "NN" in row["rna_seq_curated"] and "*[" not in row["rna_seq_curated"])
    check("每列记录目标字母表",
          all(row.get(f"{col}_alphabet") for col in candidate_names))
    check("每行 JSON 记录映射明细",
          json.loads(row["rna_seq_curation_summary"]).get("mapping"))
    check("对齐注释列同步还原且等长",
          len(row["structure_curated"]) == len(row["rna_seq_curated"]))
    check("对齐列不做字母表映射（仍保留结构符号）",
          "(" in row["structure_curated"])

    # 失败保护
    explicit_dir = os.path.join(tmp, "explicit")
    try:
        run(_parser().parse_args([
            source, "--seq-col", "Sequence,dna_seq", "--out", explicit_dir]))
        explicit_error = ""
    except CurateInputError as exc:
        explicit_error = str(exc)
    check("显式不存在列时明确失败",
          "指定列不存在" in explicit_error and "Sequence" in explicit_error)
    check("显式列错误时不生成结果",
          not os.path.exists(os.path.join(explicit_dir, "cleaned.csv")))

    single = run(_parser().parse_args([
        source, "--seq-col", "dna_seq", "--out", os.path.join(tmp, "single")]))
    check("显式单列模式不处理其他列", single["sequence_columns"] == ["dna_seq"])

    no_seq = os.path.join(tmp, "no_sequence.csv")
    with open(no_seq, "w", encoding="utf-8", newline="") as dst:
        writer = csv.DictWriter(dst, fieldnames=["id", "note"])
        writer.writeheader()
        writer.writerows({"id": r["id"], "note": r["note"]} for r in sample)
    failed_dir = os.path.join(tmp, "failed")
    stdout = io.StringIO()
    with contextlib.redirect_stdout(stdout):
        failed_code = main([no_seq, "--out", failed_dir])
    check("自动检测失败明确提示而非猜测",
          failed_code == 2 and "未能可靠识别" in stdout.getvalue())
    check("自动检测失败不生成结果",
          not os.path.exists(os.path.join(failed_dir, "cleaned.csv")))

print()
if FAILURES:
    raise SystemExit(f"失败 {len(FAILURES)} 项: {FAILURES}")
print("所有断言通过")
