"""semantic-seq-curator-agent — 生物序列语义感知清理与标注 Agent。

与 seismic cleanfa / SNIP / ShortRead / ViralEntropR 等"哑清洗"工具的本质区别：
本包不把非标准字符当作需要删除/替换的噪声，而是把每个非标准字符当作
信息载体，逐位判断其生物学语义（测序错误 / IUPAC 模糊 / 功能修饰 / 未知异常），
记录可审计的决策理由，并根据下游训练任务自动标注。

架构（真正的 Agent）：
- tools.py    Tool 抽象 + TOOL_REGISTRY（5 个可注册工具）+ 工具函数
- agent.py    Reasoner 抽象（RuleReasoner / LLMReasoner）+ 显式 ReAct 循环
- iupac.py    知识源：IUPAC 映射表、功能修饰表、保守基序库
- downstream.py 下游任务 schema
"""

from .agent import (
    SemanticSeqCuratorAgent,
    Reasoner,
    RuleReasoner,
    LLMReasoner,
    CurationRecord,
    SEMANTIC_CLASSES,
)
from .autonomous import AutonomousCurator
from .tools import (
    Tool,
    TOOL_REGISTRY,
    lookup_iupac_tool,
    check_context_tool,
    fetch_downstream_schema_tool,
    quality_evidence_tool,
    modification_knowledge_tool,
)

__all__ = [
    "SemanticSeqCuratorAgent",
    "Reasoner",
    "RuleReasoner",
    "LLMReasoner",
    "CurationRecord",
    "SEMANTIC_CLASSES",
    "AutonomousCurator",
    "Tool",
    "TOOL_REGISTRY",
    "lookup_iupac_tool",
    "check_context_tool",
    "fetch_downstream_schema_tool",
    "quality_evidence_tool",
    "modification_knowledge_tool",
]
__version__ = "0.2.0"
