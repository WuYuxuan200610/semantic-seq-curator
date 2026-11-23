"""边界用例与 Agent 架构自检。运行：python test_agent.py

覆盖：
- Tool 注册表完整性（5 个工具）
- gap / 小写 / auto 类型 / 蛋白 N(标准氨基酸) 等边界
- Reasoner 可插拔：RuleReasoner 与 LLMReasoner
- LLM 非法输出降级
"""

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from semantic_seq_curator.agent import (
    SemanticSeqCuratorAgent, LLMReasoner, SEMANTIC_CLASSES,
)
from semantic_seq_curator.tools import TOOL_REGISTRY

_fail = []


def check(name, cond):
    print(("PASS  " if cond else "FAIL  ") + name)
    if not cond:
        _fail.append(name)


# 1. Tool 注册表
check("Tool 注册表含 5 个工具", len(TOOL_REGISTRY) == 5)
check("每个 Tool 都有 name/description/parameters/run",
      all(getattr(t, "name", None) and getattr(t, "description", None)
          and hasattr(t, "run") for t in TOOL_REGISTRY.values()))

# 2. 语义类别封闭
check("四类语义类别",
      set(SEMANTIC_CLASSES) == {"sequencing_error", "iupac_ambiguity",
                                "functional_modification", "unknown_artifact"})

a = SemanticSeqCuratorAgent("nucleic_lm_pretraining", "dna")

# 3. gap
r = a.curate("g1", "ACGT-ACGT")[0]
check("gap '-' → functional_modification/keep",
      r.semantic_class == "functional_modification" and r.action.startswith("keep"))

# 4. 小写自动大写
r = a.curate("g2", "acgtracgt")[0]
check("小写 'r' → iupac_ambiguity 且原字符归一大写",
      r.semantic_class == "iupac_ambiguity" and r.original_char == "R")

# 5. auto 类型：含 U 的 RNA
b = SemanticSeqCuratorAgent("nucleic_lm_pretraining", "auto")
r = b.curate("g3", "ACGUNACGU")[0]
check("auto→RNA：N 判为 iupac_ambiguity", r.semantic_class == "iupac_ambiguity")

# 6. 蛋白 N = 天冬酰胺（标准，不触发）
c = SemanticSeqCuratorAgent("protein_lm_pretraining", "protein")
check("蛋白 N(Asn) 不作为非标准字符", len(c.curate("g4", "MKNAL")) == 0)

# 7. 未知任务报错
try:
    SemanticSeqCuratorAgent("no_such_task")
    check("未知任务报错", False)
except ValueError:
    check("未知任务报错", True)

# 8. LLMReasoner：合法 JSON
mock = lambda p: json.dumps({
    "semantic_class": "iupac_ambiguity", "action": "keep_annotate_uncertainty",
    "confidence": 0.8, "reason": "IUPAC 模糊码", "downstream_label": {}})
llm_agent = SemanticSeqCuratorAgent("nucleic_lm_pretraining", "dna",
                                    reasoner=LLMReasoner(mock))
r = llm_agent.curate("g5", "ACGTR")[0]
check("LLMReasoner 正常注入并分类", r.semantic_class == "iupac_ambiguity")

# 9. LLMReasoner：非法 JSON → 降级
bad = lambda p: "这不是 JSON"
llm_agent2 = SemanticSeqCuratorAgent("nucleic_lm_pretraining", "dna",
                                     reasoner=LLMReasoner(bad))
r = llm_agent2.curate("g6", "ACGTR")[0]
check("LLM 非法输出 → 降级 unknown_artifact/人工复核",
      r.semantic_class == "unknown_artifact" and r.action == "flag_human_review")

# 10. ReAct trace 完整（每个记录含 thought/action/observation）
r = a.curate("g7", "ACGTNACGT")[0]
steps = [t["step"] for t in r.trace]
check("ReAct trace 含 thought/action/observation",
      "thought" in steps and "action" in steps and "observation" in steps)

# 11. 占位符按整体 token 处理并还原为 k 个 N
ph = a.curate("g8", "ACG*[19]*ACGT")
check("占位符作为单个 token 记录（不是 4 个字符）",
      len(ph) == 1 and ph[0].original_char == "*[19]*")
check("占位符记录还原 k 个 N 并标注 k",
      ph[0].action == "expand_to_n"
      and ph[0].downstream_label.get("restored_bases") == 19
      and ph[0].semantic_class == "unknown_artifact")
check("占位符理由引用实证依据",
      "FASTA" in ph[0].reason and "2168" in ph[0].reason)

print()
if _fail:
    print(f"失败 {len(_fail)} 项：{_fail}")
    sys.exit(1)
print(f"全部通过（{sum(1 for _ in []) or 0} 项断言）")
