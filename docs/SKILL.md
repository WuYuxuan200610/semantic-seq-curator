---
name: semantic-seq-curator
description: 对蛋白质/核酸序列做"语义感知"清理与标注时使用。逐位判断每个非标准字符（IUPAC 模糊码、测序错误、功能修饰、未知异常）的生物学语义，记录可审计的决策理由，并按下游训练任务（预训练/分类）自动标注。绝不把非标准字符一律当噪声删除；与 seismic cleanfa、SNIP、ShortRead、ViralEntropR 等"哑清洗"工具相对。
---

# semantic-seq-curator — 生物序列语义感知清理与标注

## 定位

现有序列清洗工具（seismic cleanfa、SNIP、ShortRead、ViralEntropR 等）几乎都在做
"哑替换"或"哑删除"：把 IUPAC 模糊码统一替换为 N，或把含模糊码的序列整条删除。
它们把"不确定性"当作需要消灭的噪声。

本 Skill 反其道而行：**把每个非标准字符当作信息载体，逐位判断其生物学语义，
记录可审计的决策理由，并根据下游训练任务自动标注。**

## 语义分类框架（逐位判断，不是整条过滤）

对每一个非标准字符，必须归入以下四类之一，并给出引用具体依据的理由：

| 类别 | 含义 | 处理 | 示例 |
|---|---|---|---|
| `sequencing_error` | 明确指向测序/格式错误（低质量 Phred 对应的碱基调用、同聚物滑动、格式乱码） | 标记，可选择移除或掩码 `<mask>` | 低质量区内的 N、乱码字符 |
| `iupac_ambiguity` | 合法 IUPAC 模糊码，代表真实生物学不确定性 | **保留**，标注不确定等级；可转换为概率分布或"不确定"通道，**不**简单替换为 N | R/Y/S/W/K/M/B/D/H/V、N、X、蛋白 B/Z/J |
| `functional_modification` | 已知对应功能修饰/特殊状态的字符 | **保留**，标注修饰类型 | 蛋白 U(硒代半胱氨酸)/O(吡咯赖氨酸)/*(终止子)、gap |
| `unknown_artifact` | 无法归类的异常 | 标记为**需人工复核**，**不自动删除** | 任何查不到的符号 |

## 逐位判断所需的外部知识（工具）

Agent 不能只靠正则匹配。必须调用或检索：

1. **`lookup_iupac_tool`** — 查询字符的 IUPAC 标准定义（合法？含义？可解释集合？）。
2. **`modification_knowledge_tool`** — 查询是否为已知功能修饰/特殊状态。
3. **`check_context_tool`** — 分析该字符附近上下文（是否在保守基序、同聚物/低复杂度区）。
4. **`quality_evidence_tool`** — 拉取该位置质量证据（Phred），区分测序错误 vs 真实不确定性。
5. **`fetch_downstream_schema_tool`** — 获取下游任务（预训练/分类）对特殊字符的容忍度与标签要求。

> 若当前环境缺少某个知识源，必须如实标注"依据不足，降级为人工复核"，不能猜测。

## 逐位判断流程（ReAct 循环）

```
遇到非标准字符
  → 思考：这个字符是什么？
  → 行动：lookup_iupac_tool + modification_knowledge_tool
  → 观察：工具返回的证据
  → 思考：结合 check_context_tool + quality_evidence_tool + fetch_downstream_schema_tool
  → 行动：生成 {semantic_class, reason, action, confidence} 标注
```

降级纪律：工具查不到 → `unknown_artifact` / 人工复核，绝不编造。

## 按下游任务自动标注（标签 schema）

- **语言模型预训练**（`protein_lm_pretraining` / `nucleic_lm_pretraining`）：
  每位置输出 `{char, semantic_class, confidence, action_taken, uncertainty_level, encoding_hint}`；
  每序列输出质量标签与可用的掩码建议。
- **分类任务**（如 `enzyme_classification`）：
  输出序列级 `{usable, quality_score, mask_suggestion, functional_flag, review_required}`。

标签必须与逐位决策记录绑定，不能脱离理由单独存在。

## 决策记录格式（可审计核心）

输出结构化 `curation_log`（JSON Lines），每行对应一个非标准字符：

```json
{"seq_id":"seq_dna","position":7,"original_char":"N",
 "semantic_class":"iupac_ambiguity",
 "reason":"IUPAC 模糊码：any base (A/C/G/T)。位置 7 位于保守基序 GATA-family TF binding site 内…",
 "action":"keep_annotate_uncertainty","confidence":0.85,
 "downstream_label":{"uncertainty_level":"defined","ambiguity_set":["A","C","G","T"]}}
```

`reason` 必须引用具体依据（如"IUPAC 标准：R = A/G""位置 42 位于已知活性位点 motif 内"），
不能写"根据规则处理"。

## 依赖与实现

本 Skill 的可执行核心是纯 Python 包（无第三方依赖）：

- `semantic_seq_curator/iupac.py` — IUPAC 标准映射表、功能修饰表、保守基序库（知识源）
- `semantic_seq_curator/downstream.py` — 下游任务 schema
- `semantic_seq_curator/tools.py` — 5 个工具
- `semantic_seq_curator/agent.py` — ReAct 循环 + 语义分类策略 + `CurationRecord`
- `run_demo.py` — 最小验证（同一 N 不同上下文不同决策）

运行：

```bash
python run_demo.py   # 产出 validation/curation_log.jsonl
```

## 竞品调研与空白

详见 [RESEARCH.md](RESEARCH.md)。核心结论：现有工具是"字符级规则清洗"（哑替换/整条删除），
本 Skill 是"语义级判断与记录"；bio-rs 是唯一做逐位符号级校验的工具但无语义归因，
data-labelling 是可审计标注架构的最佳借鉴对象但对象是元数据而非序列。
