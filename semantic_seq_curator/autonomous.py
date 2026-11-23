"""AutonomousCurator —— 真正的 Agent：LLM 直接面对原始序列，自主发现语义。

与 RuleReasoner 的本质区别：
  RuleReasoner        : 规则预扫描(小写/占位符/gap 都是硬编码) → LLM 只做四分类
  AutonomousCurator   : LLM 拿到【原始序列】→ 自主观察异常 → 自主判断语义 → 决策+理由

关键：小写=低置信度、*[N]*=占位符、N=不确定、gap=对齐 这些语义，
【不是】查表查出来的，【是】LLM 用领域知识现场推断出来的。因此它能处理
规则没写过的新模式，并给出可解释的理由。

代价（必须诚实面对）：
  1. 不稳定 —— 同一输入两次推理措辞/置信度可能不同
  2. 不可复现 —— 需记录模型+提示词+参数才能审计
  3. 可能幻觉 —— 可能把 *[1]* 猜成不存在的修饰类型
  4. 成本 —— 逐条逐字符推理

正道（混合形态）：
  LLM 首次遇到新模式时自主发现语义 → 人类确认 → 沉淀为新规则 → 规则规模化执行。
  LLM 是"规则的种子来源"，规则是"LLM 的规模化零成本执行"。
"""

from __future__ import annotations

import json
from typing import Callable


# LLM 自主推理的系统提示词（核心：不喂预设规则，让 LLM 自己发现语义）
AUTONOMOUS_SYSTEM_PROMPT = (
    "你是资深生物信息学家，精通 RNA/DNA/蛋白质序列与语言模型数据整理。\n"
    "给你一段【未清洗的原始序列】，请你自主完成：\n"
    "1. 找出所有非标准字符/模式（小写字母、gap、特殊符号、模糊码、占位符等）；\n"
    "2. 对每个，判断它的生物学语义：\n"
    "   - sequencing_error    测序/格式错误（低质量、乱码）\n"
    "   - iupac_ambiguity     合法 IUPAC 模糊码（真实生物学不确定性）\n"
    "   - functional_modification 功能修饰/特殊状态（非标准氨基酸、修饰碱基、gap、终止子）\n"
    "   - low_confidence      低置信度/可变位点（如比对中的小写碱基）\n"
    "   - unknown_artifact    无法归类/已被上游破坏的占位符\n"
    "   注意：形如 *[k]* 的占位符，k 是上游丢失的碱基个数（不是编号），"
    "应还原为 k 个 N 以恢复真实长度与对齐；不要逐字符拆开处理。\n"
    "3. 决定处理动作 keep_annotate_* / mask_or_remove / expand_to_n / flag_human_review，并给理由；\n"
    "4. 铁律：绝不编造；查不到依据就标注『依据不足，人工复核』。\n"
    "只输出 JSON 数组，不要任何其它文字。\n"
)


class AutonomousCurator:
    """真 agent：把原始序列交给 LLM，LLM 自主产出 curation 记录。

    Args:
        llm_call: (prompt: str) -> str，由调用方注入（OpenAI / DSH / 任意后端）。
        downstream_task: 下游任务名（仅透传给提示词，不影响自主发现）。
    """

    def __init__(self, llm_call: Callable[[str], str], downstream_task: str = "nucleic_lm_pretraining"):
        self.llm_call = llm_call
        self.task = downstream_task
        self.learned_rules: dict[str, dict] = {}  # LLM 发现 → 人类确认后沉淀的规则

    def curate_sequence(self, seq_id: str, sequence: str, seq_type: str) -> list[dict]:
        prompt = (
            AUTONOMOUS_SYSTEM_PROMPT
            + f"\n序列类型：{seq_type}；下游任务：{self.task}。\n"
            + f"序列 id：{seq_id}；原始序列：{sequence}\n"
            + "对每个非标准字符/模式输出一条：\n"
              '{"position": <1-based 起始位置>, "original": "<字符或模式>", '
              '"semantic_class": "<五类之一>", "action": "<动作>", '
              '"reason": "<引用具体依据>", "confidence": <0~1>}'
        )
        raw = self.llm_call(prompt)
        try:
            records = json.loads(raw)
        except (json.JSONDecodeError, TypeError):
            return [{"position": None, "original": None, "semantic_class": "unknown_artifact",
                     "action": "flag_human_review", "reason": "LLM 输出非法，降级人工复核",
                     "confidence": 0.0}]
        if not isinstance(records, list):
            return [{"position": None, "original": None, "semantic_class": "unknown_artifact",
                     "action": "flag_human_review", "reason": "LLM 未返回数组", "confidence": 0.0}]
        return records

    # 人类确认后，把 LLM 发现的新语义沉淀为规则（规模化零成本执行的种子）
    def confirm_rule(self, pattern: str, semantic_class: str, action: str, note: str) -> None:
        self.learned_rules[pattern] = {
            "semantic_class": semantic_class, "action": action, "note": note,
            "source": "LLM 自主发现 + 人类确认",
        }

    def dump_learned_rules(self) -> dict:
        return dict(self.learned_rules)
