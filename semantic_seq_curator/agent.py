"""真正的 Agent：可插拔推理器 + 工具注册表 + 显式 ReAct 循环。

与"带工具注释的规则函数"的区别在于三层抽象：

1. Tool 注册表（tools.TOOL_REGISTRY）：工具是带 name/description/参数 schema 的
   一等公民，Agent 从注册表查工具、通用地调用，而不是硬编码 if-else。
2. Reasoner（推理器）：把"工具观察 → 语义决策"这一步做成可插拔接口。
   - RuleReasoner：确定性规则推理器（默认，可跑通、可审计、零依赖）
   - LLMReasoner：真正的"思考"由 LLM 完成（注入任意 llm_call 即可），
     与 AGENT.md 的 System Prompt 一一对应。
3. 显式 ReAct 循环：plan → act → observe → … → decide → emit，全程留 trace。
   这里的 plan 是确定性顺序；生产环境由 LLM 动态决定"下一步调用哪个工具"。

关键不变量：工具查不到就降级，绝不编造语义判断。
"""

from __future__ import annotations

import json
from abc import ABC, abstractmethod
from dataclasses import dataclass, field, asdict
from typing import Callable, Optional

from . import iupac as K
from . import tools as T

# 语义分类（与 System Prompt 一致）
SEMANTIC_CLASSES = (
    "sequencing_error",
    "iupac_ambiguity",
    "functional_modification",
    "unknown_artifact",
)

LOW_QUAL_THRESHOLD = 20  # Phred < 20 视为低质量

# 确定性工具计划：本实现按此顺序查证据。真正的 Agent 由 LLM 决定下一步调用哪个工具。
TOOL_PLAN = ["lookup_iupac", "modification_knowledge", "check_context", "quality_evidence"]


@dataclass
class CurationRecord:
    """curation_log 的单条记录（对应一个非标准字符）。"""
    seq_id: str
    position: int
    original_char: str
    semantic_class: str
    reason: str
    action: str
    confidence: float
    downstream_label: dict = field(default_factory=dict)  # 与逐位决策绑定的任务标注
    trace: list[dict] = field(default_factory=list)       # ReAct 追踪（审计）

    def to_jsonl(self) -> str:
        def _default(o):
            if isinstance(o, (set, frozenset)):
                return sorted(o)
            return str(o)
        return json.dumps(asdict(self), ensure_ascii=False, default=_default)


# ---------------------------------------------------------------------------
# Reasoner：观察 -> 决策 的可插拔接口
# ---------------------------------------------------------------------------
class Reasoner(ABC):
    """推理器：把工具观察 + 下游任务 schema 映射为语义决策。"""

    @abstractmethod
    def decide(self, position: int, char: str, eff_type: str,
               observations: dict, schema: dict) -> dict:
        """返回 {semantic_class, action, confidence, reason, downstream_label}。"""


class RuleReasoner(Reasoner):
    """确定性规则推理器（默认）：可跑通、可审计、零外部依赖。

    这是本 Agent 的"基准大脑"——策略公开、可复现、可单元测试。
    """

    def decide(self, position: int, char: str, eff_type: str,
               observations: dict, schema: dict) -> dict:
        return _classify(
            position, char, eff_type,
            observations["lookup_iupac"],
            observations["modification_knowledge"],
            observations["check_context"],
            observations["quality_evidence"],
            schema,
        )


class LLMReasoner(Reasoner):
    """LLM 推理器：真正的"思考-判断"由大模型完成。

    `llm_call(prompt: str) -> str` 由调用方注入（OpenAI / DSH / 任意后端）。
    输出必须是 JSON。若 LLM 返回非法 JSON，降级为人工复核，绝不猜测。
    提示词与 AGENT.md 的 System Prompt 一致。
    """

    def __init__(self, llm_call: Callable[[str], str]):
        self.llm_call = llm_call

    def decide(self, position: int, char: str, eff_type: str,
               observations: dict, schema: dict) -> dict:
        prompt = self._render_prompt(position, char, eff_type, observations, schema)
        raw = self.llm_call(prompt)
        try:
            decision = json.loads(raw)
        except (json.JSONDecodeError, TypeError):
            return {
                "semantic_class": "unknown_artifact",
                "action": "flag_human_review",
                "confidence": 0.0,
                "reason": "LLM 推理器输出非法 JSON，按纪律降级人工复核。",
                "downstream_label": {},
            }
        # 兜底：LLM 返回字段缺失时补齐
        decision.setdefault("semantic_class", "unknown_artifact")
        decision.setdefault("action", "flag_human_review")
        decision.setdefault("confidence", 0.0)
        decision.setdefault("reason", "LLM 未给出充分理由，降级。")
        decision.setdefault("downstream_label", {})
        if decision["semantic_class"] not in SEMANTIC_CLASSES:
            decision["semantic_class"] = "unknown_artifact"
        return decision

    @staticmethod
    def _render_prompt(position: int, char: str, eff_type: str,
                       observations: dict, schema: dict) -> str:
        obs = json.dumps(observations, ensure_ascii=False, default=sorted)
        return (
            "你是精通生物信息学与语言模型训练的数据管家。对序列中的非标准字符逐位判断"
            "其生物学语义，绝不默认删除。四类：sequencing_error / iupac_ambiguity / "
            "functional_modification / unknown_artifact。\n"
            f"当前字符：位置 {position}，'{char}'，序列类型 {eff_type}。\n"
            f"工具观察：{obs}\n"
            f"下游任务 schema：{json.dumps(schema, ensure_ascii=False)}\n"
            "工具查不到就降级为 unknown_artifact / 人工复核，不猜测。\n"
            "只输出 JSON：{\"semantic_class\", \"action\", \"confidence\", "
            "\"reason\", \"downstream_label\"}，reason 必须引用具体依据。"
        )


# ---------------------------------------------------------------------------
# Agent：工具注册表驱动的 ReAct 循环
# ---------------------------------------------------------------------------
class SemanticSeqCuratorAgent:
    """语义感知序列清理 Agent。

    Args:
        downstream_task: 下游任务名（见 downstream.DOWNSTREAM_SCHEMAS）。
        seq_type:       dna / rna / protein / auto（auto 自动预判）。
        quality_track:  可选 {position: phred}，用于区分测序错误与真实不确定性。
        verbose:        是否打印 ReAct 追踪。
        reasoner:       推理器（默认 RuleReasoner；可注入 LLMReasoner）。
    """

    def __init__(self, downstream_task: str, seq_type: str = "auto",
                 quality_track: Optional[dict] = None, verbose: bool = False,
                 reasoner: Optional[Reasoner] = None):
        self.task = downstream_task.strip().lower()
        self.seq_type = K.normalize_seq_type(seq_type)
        self.quality_track = quality_track or {}
        self.verbose = verbose
        self.reasoner = reasoner or RuleReasoner()
        self.downstream_schema = T.fetch_downstream_schema_tool(self.task)
        if not self.downstream_schema["found"]:
            from .downstream import known_tasks
            raise ValueError(f"未知下游任务 '{downstream_task}'，可选：{known_tasks()}")

    # ------------------------------------------------------------------
    def curate(self, seq_id: str, sequence: str) -> list[CurationRecord]:
        """对整条序列扫描，对每个非标准字符/占位符跑 ReAct 循环。

        占位符 *[k]* 必须按整体 token 处理（不是逐字符），否则 *、[、数字、]
        会被当成四个互不相关的未知字符，k 的计数含义会丢失。
        """
        import re

        raw = sequence
        seq = sequence.upper()
        eff_type = self.seq_type
        if eff_type == "auto":
            eff_type = K.infer_seq_type(seq)
        std = K.get_standard_set(eff_type)

        records: list[CurationRecord] = []
        i = 0
        while i < len(seq):
            m = re.match(K.PLACEHOLDER_PATTERN, seq[i:])
            if m:
                k = int(m.group(1))
                records.append(self._placeholder_record(
                    seq_id, i + 1, m.group(0), k, raw[i:i + m.end()]
                ))
                i += m.end()
                continue
            ch = seq[i]
            if ch not in std:
                records.append(self._react(seq_id, i + 1, ch, seq, eff_type))
            i += 1
        return records

    def _placeholder_record(self, seq_id: str, position: int, token: str,
                            k: int, original_token: str) -> CurationRecord:
        """占位符 *[k]* 的记录：原始身份不可知，但长度可还原 → k 个 N。"""
        reason = (
            f"上游占位符 {token}：原始碱基身份不可恢复，但 k={k} 是该处丢失的碱基数，"
            f"已实证（用未损坏 FASTA 侧翼锚定量真实跨度，2168 例 100% 等于 k；"
            f"标记宽度 {len(original_token)} 与丢失长度无关）。"
            f"按 k 个 N 还原：N 正确表达『未知碱基』，且长度回到真实生物学长度，"
            f"可保持与结构注释的逐位对齐并直接用于训练。"
        )
        trace = [
            {"step": "thought", "text": f"位置 {position} 发现占位符 {token}，"
                                        f"需按整体 token 处理，不能逐字符拆解。"},
            {"step": "action", "tool": "iupac.PLACEHOLDER_PATTERN",
             "args": {"token": token}},
            {"step": "observation", "result": {"k": k, "token_width": len(original_token),
                                               "semantic": "upstream_lost_bases"}},
            {"step": "thought", "text": reason},
            {"step": "action", "tool": "emit_curation_record",
             "args": {"semantic_class": "unknown_artifact", "action": "expand_to_n",
                      "confidence": 0.95}},
        ]
        return CurationRecord(
            seq_id=seq_id, position=position, original_char=token,
            semantic_class="unknown_artifact", action="expand_to_n",
            reason=reason, confidence=0.95,
            downstream_label={"restored_bases": k, "restored_as": "N",
                              "origin": "upstream_placeholder",
                              "original_identity_recoverable": False},
            trace=trace,
        )

    # ------------------------------------------------------------------
    def _react(self, seq_id: str, position: int, char: str,
               sequence: str, eff_type: str) -> CurationRecord:
        """单个非标准字符的 ReAct 循环：plan → act → observe → … → decide → emit。"""
        trace: list[dict] = []
        observations: dict = {}

        def log(kind: str, payload: dict):
            trace.append({"step": kind, **payload})
            if self.verbose:
                print(f"    [{kind}] {payload.get('text', payload)}")

        # 行动-观察：按计划从注册表查工具并调用
        for tool_name in TOOL_PLAN:
            tool = T.TOOL_REGISTRY[tool_name]
            args = self._args_for(tool_name, seq_id, position, char, sequence, eff_type)
            log("thought", {"text": f"需要证据：{tool.description}"})
            log("action", {"tool": tool.name, "args": args})
            observations[tool_name] = tool.run(**args)
            log("observation", {"result": observations[tool_name]})

        # 思考 + 决策：把观察交给推理器
        schema = self.downstream_schema["schema"]
        decision = self.reasoner.decide(position, char, eff_type, observations, schema)
        log("thought", {"text": decision["reason"]})
        log("action", {"tool": "emit_curation_record", "args": {
            "semantic_class": decision["semantic_class"],
            "action": decision["action"],
            "confidence": decision["confidence"]}})

        return CurationRecord(
            seq_id=seq_id,
            position=position,
            original_char=char,
            semantic_class=decision["semantic_class"],
            reason=decision["reason"],
            action=decision["action"],
            confidence=decision["confidence"],
            downstream_label=decision["downstream_label"],
            trace=trace,
        )

    def _args_for(self, tool_name: str, seq_id: str, position: int,
                  char: str, sequence: str, eff_type: str) -> dict:
        """为每个工具组装实参（确定性映射；真 Agent 由 LLM 产出实参）。"""
        if tool_name == "lookup_iupac":
            return {"char": char, "seq_type": eff_type}
        if tool_name == "modification_knowledge":
            return {"char": char, "seq_type": eff_type}
        if tool_name == "check_context":
            return {"seq_id": seq_id, "position": position,
                    "sequence": sequence, "seq_type": eff_type}
        if tool_name == "quality_evidence":
            return {"seq_id": seq_id, "position": position,
                    "quality_track": self.quality_track}
        if tool_name == "fetch_downstream_schema":
            return {"task": self.task}
        raise KeyError(f"未知工具 {tool_name}")


# ---------------------------------------------------------------------------
# 规则推理器的决策策略（模块级，可单测）
# ---------------------------------------------------------------------------
def _decide(semantic_class, action, confidence, reason, downstream_label) -> dict:
    return {
        "semantic_class": semantic_class,
        "action": action,
        "confidence": confidence,
        "reason": reason,
        "downstream_label": downstream_label,
    }


def _classify(position, char, eff_type, lookup, mod, ctx, qual, schema) -> dict:
    """策略层：把工具观察映射为语义分类 + 动作 + 置信度 + 理由。"""

    # (1) 功能修饰 / 特殊状态
    if mod["found"]:
        mtype = mod["modification_type"]
        reason = (f"功能修饰/特殊状态：{mod['meaning']}。依据 IUPAC-JCBN 标准，"
                  f"'{char}' 对应 {mtype}，属真实生物学信息，保留并标注修饰类型。")
        return _decide("functional_modification", "keep_annotate_modification",
                       0.95, reason, {"modification_type": mtype})

    # (2) 比对 gap
    if lookup.get("is_gap"):
        reason = (f"比对 gap（'{char}'）：代表已知缺失/位置状态，而非身份不确定性，"
                  f"保留并标注为 gap。")
        return _decide("functional_modification", "keep_annotate_gap",
                       0.9, reason, {"modification_type": "gap"})

    # (3) 合法 IUPAC 模糊码
    if lookup["found"] and lookup["is_ambiguity"]:
        low_quality = qual["found"] and qual["phred"] < LOW_QUAL_THRESHOLD
        if low_quality and ctx["homopolymer"]:
            reason = (f"测序错误：'{char}' 虽是合法 IUPAC 码，但质量证据 Phred={qual['phred']} "
                      f"(<{LOW_QUAL_THRESHOLD}) 且处于同聚物/低复杂度区 "
                      f"(侧翼 {ctx['flanking']})，更可能为低置信调用/同聚物滑动，建议掩码或移除。")
            return _decide("sequencing_error", "mask_or_remove", 0.7, reason,
                           {"uncertainty_level": "high"})
        # 真实生物学不确定性
        if ctx["in_motif"]:
            motif_names = "、".join(m["name"] for m in ctx["motifs"])
            reason = (f"IUPAC 模糊码：{lookup['meaning']}。位置 {position} 位于保守基序 "
                      f"{motif_names} 内（侧翼 {ctx['flanking']}），该不确定性在功能位点上有意义，"
                      f"保留并标注不确定等级。")
            return _decide("iupac_ambiguity", "keep_annotate_uncertainty", 0.85, reason,
                           _ambig_label(lookup))
        conf = 0.6 if not qual["found"] else 0.75
        quality_note = (f"质量证据 Phred={qual['phred']}。" if qual["found"]
                        else "无质量证据，依据不足但码本身合法。")
        reason = (f"IUPAC 模糊码：{lookup['meaning']}。{quality_note}"
                  f"代表真实生物学不确定性，保留并标注不确定等级；"
                  f"可转换为概率分布而非简单替换为 N。")
        return _decide("iupac_ambiguity", "keep_annotate_uncertainty", conf, reason,
                       _ambig_label(lookup))

    # (4) 无法归类的异常
    if qual["found"] and qual["phred"] < LOW_QUAL_THRESHOLD:
        reason = (f"测序错误（存疑）：'{char}' 未在 IUPAC 表中查到，"
                  f"但质量证据 Phred={qual['phred']} 极低，疑似格式乱码/低置信调用，"
                  f"建议掩码并保留日志供复核。")
        return _decide("sequencing_error", "mask_or_remove", 0.5, reason,
                       {"uncertainty_level": "max"})
    reason = (f"依据不足：'{char}' 未在 IUPAC 标准表或修饰知识源中查到，"
              f"且无质量证据支撑任何解释。按纪律降级为人工复核，不自动删除、不猜测。")
    return _decide("unknown_artifact", "flag_human_review", 0.3, reason, {})


def _ambig_label(lookup: dict) -> dict:
    return {
        "uncertainty_level": "defined" if lookup["canonical_set"] else "max",
        "ambiguity_set": sorted(lookup["canonical_set"]) if lookup["canonical_set"] else None,
    }
