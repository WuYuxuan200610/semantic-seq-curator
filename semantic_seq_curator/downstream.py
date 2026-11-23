"""下游训练任务 schema。

不同下游任务对"非标准字符"的容忍度与标签要求不同，Agent 的逐位决策
（保留/掩码/删除/转换）必须随任务变化。本模块是 fetch_downstream_schema_tool
返回内容的知识源。
"""

from __future__ import annotations

DOWNSTREAM_SCHEMAS: dict[str, dict] = {
    "protein_lm_pretraining": {
        "task": "蛋白质语言模型预训练（如 ESM-2 / ProtBERT 风格）",
        "tolerance_policy": (
            "iupac_ambiguity(B/Z/J/X) 保留为独立「不确定」通道或 one-hot 概率分布，"
            "不折叠成单一残基；functional_modification(U/O/*) 保留为独立 token（U/O）"
            "或按 terminator 处理(*)；sequencing_error 建议掩码后按 mask token 训练；"
            "unknown_artifact 一律不进入训练集，转人工复核。"
        ),
        "label_schema": {
            "per_position": ["char", "semantic_class", "confidence", "action_taken",
                             "uncertainty_level", "encoding_hint"],
            "per_sequence": ["seq_id", "usable", "n_nonstandard", "n_error", "n_ambiguity",
                             "n_modification", "n_unknown", "quality_score"],
        },
        "mask_suggestion": "对 sequencing_error 位置输出 <mask>；对 X 输出 <unk> 通道。",
    },
    "nucleic_lm_pretraining": {
        "task": "核酸语言模型预训练（如 DNABERT / Nucleotide Transformer 风格）",
        "tolerance_policy": (
            "合法 IUPAC 模糊码(R/Y/S/W/K/M/B/D/H/V/N) 保留并记录不确定集合，"
            "可编码为概率分布而非单一 token；X 保留为 <unk> 通道；"
            "修饰碱基按约定 token 保留；sequencing_error 掩码；unknown 转人工。"
        ),
        "label_schema": {
            "per_position": ["char", "semantic_class", "confidence", "action_taken",
                             "ambiguity_set", "encoding_hint"],
            "per_sequence": ["seq_id", "usable", "n_nonstandard", "n_error", "n_ambiguity",
                             "n_modification", "n_unknown", "quality_score"],
        },
        "mask_suggestion": "对 sequencing_error 输出 <mask>；对 N/X 输出概率向量通道。",
    },
    "enzyme_classification": {
        "task": "酶功能分类（监督式序列分类，如 EC number 预测）",
        "tolerance_policy": (
            "分类任务对单个模糊位置较敏感：保守区内的 iupac_ambiguity 保留并提示"
            "分类器用不确定通道；非保守区低置信度的 sequencing_error 建议掩码或移除；"
            "functional_modification(U/O) 保留为真实残基信号；unknown 转人工。"
        ),
        "label_schema": {
            "per_sequence": ["seq_id", "usable", "quality_score", "mask_suggestion",
                             "functional_flag", "review_required"],
        },
        "mask_suggestion": "低置信 sequencing_error 位置掩码；保守区模糊位置保留并标注。",
    },
}


def fetch_downstream(task: str) -> dict:
    """按任务名取 schema；未知任务返回 None，Agent 据此降级。"""
    return DOWNSTREAM_SCHEMAS.get((task or "").strip().lower())


def known_tasks() -> list[str]:
    return list(DOWNSTREAM_SCHEMAS.keys())
