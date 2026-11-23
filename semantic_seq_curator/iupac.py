"""IUPAC 标准映射表与外部知识源（知识源，非正则硬编码）。

Agent 的逐位判断依赖这些"外部知识"，而非简单的正则匹配。当某个字符
在下列知识源中都查不到时，Agent 必须降级为 `unknown_artifact`（人工复核），
绝不编造语义。

来源约定（供 reason 字段引用）：
- 核酸 IUPAC 码：IUPAC-IUB 标准 (Cornish-Bowden 1985, NAR 13:3021)。
- 蛋白质 IUPAC 码：IUPAC-IUB Joint Commission on Biochemical Nomenclature。
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# 核酸 IUPAC 码（DNA/RNA）
# ---------------------------------------------------------------------------
DNA_STANDARD = frozenset("ACGT")
RNA_STANDARD = frozenset("ACGU")

# set 表示该码可被解释为的碱基集合；None 表示"未知/掩蔽"，无确定集合。
NUC_IUPAC: dict[str, dict] = {
    "A": {"meaning": "Adenine", "set": frozenset("A"), "standard": True},
    "C": {"meaning": "Cytosine", "set": frozenset("C"), "standard": True},
    "G": {"meaning": "Guanine", "set": frozenset("G"), "standard": True},
    "T": {"meaning": "Thymine", "set": frozenset("T"), "standard": True},
    "U": {"meaning": "Uracil", "set": frozenset("U"), "standard": True},
    "R": {"meaning": "purine (A/G)", "set": frozenset("AG"), "ambiguity": True},
    "Y": {"meaning": "pyrimidine (C/T)", "set": frozenset("CT"), "ambiguity": True},
    "S": {"meaning": "strong (G/C)", "set": frozenset("GC"), "ambiguity": True},
    "W": {"meaning": "weak (A/T)", "set": frozenset("AT"), "ambiguity": True},
    "K": {"meaning": "keto (G/T)", "set": frozenset("GT"), "ambiguity": True},
    "M": {"meaning": "amino (A/C)", "set": frozenset("AC"), "ambiguity": True},
    "B": {"meaning": "not A (C/G/T)", "set": frozenset("CGT"), "ambiguity": True},
    "D": {"meaning": "not C (A/G/T)", "set": frozenset("AGT"), "ambiguity": True},
    "H": {"meaning": "not G (A/C/T)", "set": frozenset("ACT"), "ambiguity": True},
    "V": {"meaning": "not T (A/C/G)", "set": frozenset("ACG"), "ambiguity": True},
    "N": {"meaning": "any base (A/C/G/T)", "set": frozenset("ACGT"), "ambiguity": True},
    # X：非严格 IUPAC 标准码，但数据库/比对中普遍用作"掩蔽/未知"。
    "X": {"meaning": "masked/unknown base", "set": None, "ambiguity": True,
          "note": "non-standard IUPAC; commonly 'masked/unknown'"},
    "-": {"meaning": "alignment gap", "set": None, "is_gap": True},
    ".": {"meaning": "alignment gap", "set": None, "is_gap": True},
}

# ---------------------------------------------------------------------------
# 蛋白质 IUPAC 码
# ---------------------------------------------------------------------------
PROTEIN_STANDARD = frozenset("ARNDCQEGHILKMFPSTWYV")

PROTEIN_IUPAC: dict[str, dict] = {
    "B": {"meaning": "Asx (D/N)", "set": frozenset("DN"), "ambiguity": True},
    "Z": {"meaning": "Glx (E/Q)", "set": frozenset("EQ"), "ambiguity": True},
    "J": {"meaning": "Xle (I/L)", "set": frozenset("IL"), "ambiguity": True},
    "X": {"meaning": "unknown/any residue", "set": None, "ambiguity": True},
    "U": {"meaning": "Selenocysteine (Sec), 21st amino acid", "modification": "selenocysteine"},
    "O": {"meaning": "Pyrrolysine (Pyl), 22nd amino acid", "modification": "pyrrolysine"},
    "*": {"meaning": "translation stop / terminator", "modification": "stop"},
    "-": {"meaning": "alignment gap", "is_gap": True},
    ".": {"meaning": "alignment gap", "is_gap": True},
}

# ---------------------------------------------------------------------------
# 功能修饰 / 特殊状态（补充知识源，供 modification_knowledge_tool 查询）
# 核酸的修饰碱基（m6A、5mC、ψ 等）在裸 FASTA 中通常不以单字母出现；
# 若下游数据用约定字母表示，则登记于此，未登记的降级为 unknown_artifact。
# ---------------------------------------------------------------------------
NUC_MODIFICATIONS: dict[str, dict] = {
    # 示例约定（如出现则保留并标注）；默认不登记未知符号。
}

# ---------------------------------------------------------------------------
# 保守基序库（上下文判断的本地知识源，非穷举；查不到即"依据不足"）
# ---------------------------------------------------------------------------
MOTIF_DB: dict[str, list[dict]] = {
    "dna": [
        {"name": "TATA box (core promoter)", "pattern": "TATAAA",
         "role": "transcription initiation"},
        {"name": "polyadenylation signal", "pattern": "AATAAA",
         "role": "mRNA 3'-end processing"},
        {"name": "GATA-family TF binding site", "pattern": "GATAA",
         "role": "transcription factor binding (conserved)"},
        {"name": "Kozak-like context", "pattern": "GCCRCC",
         "role": "translation initiation (R=A/G)"},
    ],
    "rna": [
        {"name": "Shine-Dalgarno", "pattern": "AGGAGG",
         "role": "ribosome binding"},
        {"name": "branch point consensus", "pattern": "YNCURAY",
         "role": "splicing branch site"},
    ],
    "protein": [
        {"name": "P-loop (Walker A)", "pattern": "GXXXXGK[TS]",
         "role": "ATP/GTP nucleotide binding"},
        {"name": "Ser-His-Asp catalytic triad context", "pattern": "GXSXG",
         "role": "hydrolase active site"},
        {"name": "DEAD-box helicase motif", "pattern": "DEAD",
         "role": "RNA helicase"},
        {"name": "Cys-x-x-Cys (redox)", "pattern": "CXXC",
         "role": "disulfide/redox active site"},
    ],
}

# ---------------------------------------------------------------------------
# 上游占位符：*[k]*（k = 该处被替换掉的原始碱基个数）
#
# 语义已实证（非假设）：用未损坏的 FASTA_sequence 对占位符两侧做锚定，
# 量出两锚点之间的真实跨度，2168 例中跨度等于 k 的比例为 100%
# （锚点≥15bp；10bp 时因误匹配有 3/2967 偏差）。因此：
#   *[k]*  →  k 个 N   （还原真实长度，恢复与结构注释的逐位对齐）
# 注意：标记文本自身宽度（如 *[19]* 宽 6）与 k 无关，不能用它当长度；
# Structure 列里的 '~' 个数等于标记宽度，是从损坏序列算出的循环证据。
# ---------------------------------------------------------------------------
PLACEHOLDER_PATTERN = r"\*\[ *(\d+)\]\*"
PLACEHOLDER_NOTE = (
    "*[k]* 表示上游丢失了 k 个碱基，应还原为 k 个 N；"
    "k 已经 FASTA 侧翼锚定实证，标记宽度不代表缺失长度。"
)


def normalize_seq_type(seq_type: str) -> str:
    """把用户给的序列类型归一化为 dna / rna / protein / auto。"""
    t = (seq_type or "auto").strip().lower()
    mapping = {
        "dna": "dna", "rna": "rna", "protein": "protein", "prot": "protein",
        "aa": "protein", "peptide": "protein", "auto": "auto",
    }
    return mapping.get(t, "auto")


def infer_seq_type(sequence: str) -> str:
    """基于整条序列的字母表归属做粗粒度类型预判（借 dsh-bio-genie 思路，仅作前置）。"""
    seq = sequence.upper()
    letters = set(seq)
    if "U" in letters and "T" not in letters:
        return "rna"
    if letters & set("ACGTUNRYSWKMBDHVX-."):
        # 出现核酸特征字母但无非核酸字母 -> 核酸
        if letters <= set("ACGTUNRYSWKMBDHVX-."):
            return "rna" if ("U" in letters) else "dna"
    return "protein"


def get_iupac_table(seq_type: str) -> dict[str, dict]:
    """返回该序列类型对应的 IUPAC 表。"""
    t = normalize_seq_type(seq_type)
    if t == "protein":
        return PROTEIN_IUPAC
    return NUC_IUPAC


def get_standard_set(seq_type: str) -> frozenset:
    t = normalize_seq_type(seq_type)
    if t == "protein":
        return PROTEIN_STANDARD
    if t == "rna":
        return RNA_STANDARD
    return DNA_STANDARD
