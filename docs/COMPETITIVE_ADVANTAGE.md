# 竞品共性 vs 本 Agent 差异化优势

> 基于 9 个工具（seismic cleanfa / SNIP / ViralEntropR / ShortRead / dsh-bio-genie /
> data-labelling / bio-rs / atomworks.io / Data-Juicer）的源码级调研，见 [RESEARCH.md](RESEARCH.md)。

## 一、竞品共性：它们"都做"与"都不做"的同一件事

把 9 个工具摊开看，虽然表面不同（替换 / 删除 / 过滤 / 校验 / 标注元数据），
但在"非标准字符"这件事上，它们共享 6 条共性：

| # | 竞品共性 | 具体表现（实例） |
|---|---|---|
| C1 | **决策是全局规则，无上下文感知** | cleanfa"所有模糊码→N"；ShortRead"含一个模糊码就删整条"；它们不看这个字符在保守区还是低质量区 |
| C2 | **把非标准字符当"噪声"，动作只有删/换/滤** | SNIP 整条删除、cleanfa 替换为 N、Data-Juicer 超比例丢弃；没有"保留并标注"这个动作 |
| C3 | **不做语义分类** | 不区分"测序错误 vs IUPAC 模糊 vs 功能修饰 vs 未知异常"；唯一例外 bio-rs 只做"合法/非法符号"校验，不归因生物学含义 |
| C4 | **不输出决策理由（不可审计）** | cleanfa/SNIP/ShortRead 只输出清洗后序列或计数；唯一例外 data-labelling 有 rationale，但对象是样本元数据 |
| C5 | **无下游任务感知** | 处理策略不随"预训练 vs 分类"变化，一套规则通吃 |
| C6 | **无降级机制，查不到就硬处理** | cleanfa 遇非法字符直接报错中止；SNIP 直接删；而不是"标记人工复核" |

**一句话共性**：它们都在做**"字符级规则清洗"**——把不确定性当噪声消灭，把决策写死成一条全局规则，
不解释、不审计、不随任务变化。

## 二、本 Agent 的差异化优势（一一对应）

| 竞品共性 | 本 Agent 的差异 |
|---|---|
| C1 全局规则、无上下文 | **逐位 + 上下文感知**：查基序、同聚物、质量分，同一个 `N` 在保守区（保留）与低质量同聚物区（掩码）给不同决策 |
| C2 当噪声、只删/换/滤 | **目标字母表驱动**：不是"尽量不删"，而是**保证输出 100% 是模型能读的 token**——能保留的按语义映射，占位符按计数还原长度；删除只是可选策略（`--gaps drop`） |
| C3 无语义分类 | **逐字符语义归因**：小写=低置信度、gap=比对缺失、`*[k]*`=上游丢失 k 个碱基、模糊码=真实不确定性，映射理由逐条记录 |
| C4 不可审计 | **可审计映射**：每条 `{原来是什么 → 映射成什么 → 为什么}`，输出经独立脚本复核合法性 |
| C5 无任务感知 | **字母表随 tokenizer 切换**：`dna`/`rna`/`protein`/`dna-iupac`/`protein-extended`，同一份数据可产出不同词表版本 |
| C6 无降级 | **降级纪律**：无法归类 → 映射为 N/X 保证可输入，并标注；绝不编造语义 |

> **诚实的边界**：在**严格字母表**下把 IUPAC 模糊码映射为 `N`，与 seismic cleanfa 一致——
> 对固定词表这是唯一正确做法，**因此"模糊码→N"本身不构成差异**。
> 真正的差异是：**① 独立复核的 100% 合法性保证；② 占位符按计数还原，把 cleanfa/哑程序丢弃的
> 序列长度补回来（数据产出 100.2% vs 93.5%）；③ 映射逐条可审计；④ 字母表可切换**。

## 三、结构性差异：真正的 Agent，不是静态流水线

竞品是"静态规则"（cleanfa/SNIP/ShortRead）或"静态 Prompt"（data-labelling 是元数据标注管道）。
本 Agent 在架构上也是 agentic 的：

| 维度 | 竞品 | 本 Agent |
|---|---|---|
| 工具 | 硬编码函数 | **Tool 注册表**（`Tool` 带 name/description/参数 schema/run，5 个工具） |
| 决策 | 写死规则 | **可插拔 Reasoner**：`RuleReasoner`（确定性、可测）+ `LLMReasoner`（真 LLM 思考，注入任意 `llm_call`） |
| 循环 | 单步变换 | **显式 ReAct 循环**：plan → act → observe → … → decide → emit，全程留 trace |
| 边界 | 报错/删除 | **降级人工复核**，非法 LLM 输出也降级 |

代码证据：
- `semantic_seq_curator/tools.py` 的 `Tool` / `TOOL_REGISTRY`（工具调用协议）
- `semantic_seq_curator/agent.py` 的 `Reasoner` / `RuleReasoner` / `LLMReasoner` + `_react` ReAct 循环
- `test_agent.py` 验证了 Tool 注册表、Reasoner 可插拔、LLM 降级、gap/小写/auto/蛋白 N 等边界

## 四、可量化的差异（金标准 11 字符）

| 方法 | 端到端正确率 | 信息保留率 | 噪声清除率 | 可解释率 |
|---|---|---|---|---|
| 白板 Agent | 18% | 11% | 50% | 100% |
| 简单 Prompt | 45% | 44% | 50% | 91% |
| seismic cleanfa | 18% | 11% | 50% | 0% |
| SNIP | 18% | 0% | 100% | 0% |
| **本 Agent** | **100%** | **100%** | **100%** | **100%** |

![对比图](../benchmark/comparison.png)
