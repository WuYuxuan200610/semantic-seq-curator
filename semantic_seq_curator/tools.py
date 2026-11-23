"""Agent 工具箱：5 个工具，每个工具返回结构化证据（含 found/sufficient 标志）。

设计原则：
- 工具只返回"查到的事实"，不做决策。决策由 agent.py 的策略层完成。
- 每个工具都必须诚实报告"查没查到"（found）与"证据是否充分"（sufficient），
  查不到时 Agent 才能据此降级为 unknown_artifact / 人工复核，绝不猜测。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

from . import iupac as K
from . import downstream as DS

# 核酸里 N/X 是"任意/未知"（可作基序匹配通配符）；蛋白质里 N 是标准天冬酰胺，
# X 才是通配符。故通配符集合随序列类型变化。
_WILDCARDS = {"dna": frozenset("NX"), "rna": frozenset("NX"), "protein": frozenset("X")}


# ---------------------------------------------------------------------------
# 工具 1：lookup_iupac_tool —— 查询字符的 IUPAC 标准定义
# ---------------------------------------------------------------------------
def lookup_iupac_tool(char: str, seq_type: str) -> dict:
    """查询字符 char 在 seq_type 字母表中的 IUPAC 标准定义。"""
    ch = char.upper()
    table = K.get_iupac_table(seq_type)
    std = K.get_standard_set(seq_type)
    if ch in table:
        entry = table[ch]
        return {
            "found": True,
            "is_standard": bool(entry.get("standard")) or ch in std,
            "is_ambiguity": bool(entry.get("ambiguity")),
            "is_gap": bool(entry.get("is_gap")),
            "meaning": entry.get("meaning", ""),
            "canonical_set": entry.get("set"),
            "note": entry.get("note", ""),
            "sufficient": True,
        }
    return {
        "found": False,
        "is_standard": False,
        "is_ambiguity": False,
        "is_gap": False,
        "meaning": f"未在 IUPAC 标准表（{K.normalize_seq_type(seq_type)}）中查到 '{char}'",
        "canonical_set": None,
        "note": "",
        "sufficient": False,
    }


# ---------------------------------------------------------------------------
# 工具 2：modification_knowledge_tool —— 查功能修饰 / 特殊状态
# ---------------------------------------------------------------------------
def modification_knowledge_tool(char: str, seq_type: str) -> dict:
    """查询字符是否为已知功能修饰 / 特殊状态（非标准氨基酸、修饰碱基、终止子等）。"""
    ch = char.upper()
    table = K.get_iupac_table(seq_type)
    entry = table.get(ch)
    if entry and entry.get("modification"):
        return {
            "found": True,
            "modification_type": entry["modification"],
            "meaning": entry.get("meaning", ""),
            "sufficient": True,
        }
    return {"found": False, "modification_type": None, "meaning": "", "sufficient": False}


# ---------------------------------------------------------------------------
# 工具 3：check_context_tool —— 分析字符附近序列上下文
# ---------------------------------------------------------------------------
def _match_at(pattern: str, window: str, seq_type: str) -> bool:
    """IUPAC 感知的基序匹配：pattern 的模糊码按其集合展开；窗口中的通配符
    （核酸 N/X、蛋白 X）可匹配 pattern 的任意字符。"""
    wc = _WILDCARDS.get(K.normalize_seq_type(seq_type), frozenset("X"))
    p, w = pattern.upper(), window.upper()
    if len(p) != len(w):
        return False
    j = 0
    while j < len(p):
        pch = p[j]
        if pch == "[":
            end = p.find("]", j)
            if end == -1:
                return False
            allowed = set(p[j + 1:end])
            if w[j] in wc or w[j] in allowed:
                j = end + 1
                continue
            return False
        if pch == "X":  # 显式通配
            j += 1
            continue
        entry = K.get_iupac_table(seq_type).get(pch)
        if entry and entry.get("set"):
            if w[j] in wc or w[j] in entry["set"]:
                j += 1
                continue
            return False
        # 普通残基：pattern 字符需与窗口字符一致，或窗口为通配符
        if w[j] in wc or w[j] == pch:
            j += 1
            continue
        return False
    return True


def _motif_spans(pattern: str, seq_window: str, seq_type: str) -> list[tuple[int, int]]:
    """返回 seq_window 内 pattern 命中的所有 [start, end) 区间（0-based，半开）。"""
    p = pattern.upper()
    s = seq_window.upper()
    spans: list[tuple[int, int]] = []
    if len(p) > len(s):
        return spans
    for i in range(len(s) - len(p) + 1):
        if _match_at(p, s[i:i + len(p)], seq_type):
            spans.append((i, i + len(p)))
    return spans


def _homopolymer_run(seq: str, pos0: int) -> int:
    """返回"包含 pos0"的同聚物游程长度（pos0 本身可能是嵌在 run 中的非标准字符）。"""
    if not (0 <= pos0 < len(seq)):
        return 0
    left_ch = seq[pos0 - 1] if pos0 > 0 else ""
    right_ch = seq[pos0 + 1] if pos0 < len(seq) - 1 else ""
    anchors = {c for c in (left_ch, right_ch) if c}
    best = 0
    for anchor in anchors:
        run = 0
        i = pos0 - 1
        while i >= 0 and seq[i] == anchor:
            run += 1
            i -= 1
        i = pos0 + 1
        while i < len(seq) and seq[i] == anchor:
            run += 1
            i += 1
        best = max(best, run)
    return best


def check_context_tool(seq_id: str, position: int, sequence: str,
                       seq_type: str, window: int = 12) -> dict:
    """分析 position（1-based）附近上下文：保守基序命中、同聚物/低复杂度信号。"""
    seq = sequence.upper()
    pos0 = position - 1
    lo = max(0, pos0 - window)
    hi = min(len(seq), pos0 + window + 1)
    flanking = seq[lo:hi]

    t = K.normalize_seq_type(seq_type)
    offset = pos0 - lo  # 目标字符在 flanking 中的 0-based 位置
    motifs = []
    for m in K.MOTIF_DB.get(t, []):
        for start, end in _motif_spans(m["pattern"], flanking, t):
            motifs.append({
                "name": m["name"],
                "role": m["role"],
                "pattern": m["pattern"],
                "covers_position": start <= offset < end,
            })

    in_motif = any(m["covers_position"] for m in motifs)
    homopolymer = _homopolymer_run(seq, pos0) >= 3

    return {
        "found": bool(motifs) or homopolymer,
        "in_motif": in_motif,
        "motifs": motifs,
        "homopolymer": homopolymer,
        "flanking": flanking,
        "sufficient": True,
    }


# ---------------------------------------------------------------------------
# 工具 4：quality_evidence_tool —— 拉取每位置质量证据（Phred）
# ---------------------------------------------------------------------------
def quality_evidence_tool(seq_id: str, position: int,
                          quality_track: dict | None = None) -> dict:
    """返回某位置的测序质量证据（Phred）。无质量数据时 found=False。"""
    if quality_track and position in quality_track:
        q = quality_track[position]
        interp = ("低质量（可能为测序错误/低置信调用）" if q < 20
                  else "高质量（模糊码更可能为真实信号）")
        return {"found": True, "phred": q, "interpretation": interp, "sufficient": True}
    return {"found": False, "phred": None, "interpretation": "无质量数据，依据不足", "sufficient": False}


# ---------------------------------------------------------------------------
# 工具 5：fetch_downstream_schema_tool —— 获取下游任务容忍度与标签要求
# ---------------------------------------------------------------------------
def fetch_downstream_schema_tool(task: str) -> dict:
    """获取下游任务对特殊字符的容忍度与标签 schema。未知任务 found=False。"""
    schema = DS.fetch_downstream(task)
    if schema:
        return {"found": True, "schema": schema, "sufficient": True}
    return {
        "found": False,
        "schema": None,
        "sufficient": False,
        "note": f"未知下游任务 '{task}'，可选：{DS.known_tasks()}",
    }


# ---------------------------------------------------------------------------
# Tool 抽象与注册表：把工具封装为"可被 Agent 动态调用"的一等公民，
# 每个工具带 name / description / 参数 schema / run。这是"真 Agent"的
# 工具调用协议（与 DSH 工具 / MCP 工具同构），Agent 从注册表查工具并调用。
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Tool:
    name: str
    description: str
    parameters: dict          # JSON-Schema 风格的参数说明（人机可读）
    run: Callable[..., dict]


TOOL_REGISTRY: dict[str, Tool] = {
    "lookup_iupac": Tool(
        name="lookup_iupac",
        description="查询字符的 IUPAC 标准定义（是否合法、含义、可解释集合）",
        parameters={"char": "string", "seq_type": "dna|rna|protein"},
        run=lookup_iupac_tool,
    ),
    "modification_knowledge": Tool(
        name="modification_knowledge",
        description="查询字符是否为已知功能修饰/特殊状态（如 U/O/*/gap）",
        parameters={"char": "string", "seq_type": "dna|rna|protein"},
        run=modification_knowledge_tool,
    ),
    "check_context": Tool(
        name="check_context",
        description="分析该字符附近上下文：是否在保守基序、是否同聚物/低复杂度区",
        parameters={"seq_id": "string", "position": "int(1-based)",
                    "sequence": "string", "seq_type": "dna|rna|protein"},
        run=check_context_tool,
    ),
    "quality_evidence": Tool(
        name="quality_evidence",
        description="拉取该位置的质量证据（Phred），区分测序错误与真实不确定性",
        parameters={"seq_id": "string", "position": "int(1-based)",
                    "quality_track": "dict|None"},
        run=quality_evidence_tool,
    ),
    "fetch_downstream_schema": Tool(
        name="fetch_downstream_schema",
        description="获取下游任务（预训练/分类）对特殊字符的容忍度与标签 schema",
        parameters={"task": "string"},
        run=fetch_downstream_schema_tool,
    ),
}
