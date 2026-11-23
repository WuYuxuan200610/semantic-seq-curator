# 竞品调研报告与空白分析

> 本报告所有结论均来自实际检索（web_search / web_fetch / 源码抓取），每条带可点击来源，未凭记忆编造。

## 一、逐工具核实

### 1. seismic cleanfa（SEISMIC-RNA 子命令）
- 仓库 [rouskinlab/seismic-rna](https://github.com/rouskinlab/seismic-rna) · 文档 [cleanfa](https://rouskinlab.github.io/seismic-rna/use/utility/cleanfa.html) · **14★**
- 核心逻辑：修整 FASTA 使其兼容 SEISMIC-RNA（去头部元数据、大写、U↔T、模糊码处理）。
- **对非标准字符**：模糊 IUPAC 码统一**替换为 N**（哑替换）；真正非字母表字符**报错中止**，要求手动删除整行。
- 逐位判断：检测层是（正则逐位替换），动作层是"统一替换/报错"；**不做生物学意义判断**；**不记录理由**；**无语义感知**。

### 2. SNIP（BBCMdP/SNIP）
- 仓库 [BBCMdP/SNIP](https://github.com/BBCMdP/SNIP) · **0★**
- 核心逻辑：`Sequence Non-IUPAC Purge script`——移除 multifasta 中含"非 IUPAC 字符"的**整条**序列，并剥离末端 `*`/`.`。
- **对非标准字符**：白名单仅 20 种标准氨基酸（`--aligned` 才额外允许 `-`/`.`），含 B/J/X/Z/U/O 等模糊码的序列被**整条删除**。
- 逐位判断：检测层逐字符找首个违规字符，动作层**整条删除**；**无"保留有意义模糊码"选项**；理由记录仅为计数（Total/Accepted/Removed），无逐残基理由；**无语义感知**。

### 3. ViralEntropR（R 包）
- CRAN [ViralEntropR](https://cran.r-project.org/web/packages/ViralEntropR/) · 仓库 [vadimtyuryaev/ViralEntropR](https://github.com/vadimtyuryaev/ViralEntropR) · 源码 [filter_ambiguous_sequences.R](https://rdrr.io/cran/ViralEntropR/src/R/filter_ambiguous_sequences.R) · **GitHub 0★**
- 核心逻辑：熵驱动病毒变异检测前处理，`filter_ambiguous_sequences()` 移除含模糊残基的序列行。
- **对非标准字符**：含至少一个 B/J/X/Z 的**整行删除**（`flag for removal`）；唯一语义区分是**保留 gap（-）**（"gap=已知缺失而非身份不确定"）。
- 逐位判断：检测层向量化逐位，动作层整行删除；理由仅返回 `DeletedSeqId / NumberAmbiguous / RangeAmbiguous`，无逐位"为什么"；**语义感知部分**（gap vs 模糊码二分类，但不区分不同模糊码）。

### 4. ShortRead（Bioconductor）
- 手册 [clean()](https://rdrr.io/bioc/ShortRead/man/clean.html) · 镜像 [Bioconductor/ShortRead](https://github.com/Bioconductor/ShortRead) · **8★**
- 核心逻辑：FASTQ 输入/过滤/质检；`clean()` 移除含模糊碱基的短读。
- **对非标准字符**：`clean` "removes all sequences containing 1 or more ambiguous bases（IUPAC symbols different from A,T,G,C）"——**整条删除**，all-or-nothing。
- 无阈值选项、无逐位替换、不记录理由、**无语义感知**。另有 `srFilter/nFilter/alphabetScore` 属数量/质量层，非语义分类。

### 5. dsh-bio-genie（DSH 插件）
- 仓库 [moonbowterfly/dsh-bio-genie](https://github.com/moonbowterfly/dsh-bio-genie) · npm [@dsh-bio/dsh-bio-genie](https://cdn.jsdelivr.net/npm/@dsh-bio/dsh-bio-genie@0.6.47/README.md) · **4★**
- 核心逻辑："许愿式"生信分析插件（57 个语义工具 + Biopython 执行器）。
- **对非标准字符**：识别 IUPAC 模糊码（R/Y/S/W/K/M/B/D/H/V）、X（未知/修饰碱基）、gap（-/.）**仅用于整条序列的类型自动判断**（含 U 无 T→RNA；含模糊码→DNA；非核酸字母→蛋白质），**不做清洗决策**。
- **无逐位语义清理能力**；理由仅针对分析结论（`bio_log`），非逐位清理理由；语义感知弱（字母表归属）。

### 6. data-labelling（nextflow-io/nextflow 示例）
- 目录 [examples/agents/14_data-labelling](https://github.com/nextflow-io/nextflow/tree/master/examples/agents/14_data-labelling) · README [data-labelling](https://github.com/nextflow-io/nextflow/blob/master/examples/agents/14_data-labelling/README.md) · 宿主 [nextflow-io/nextflow](https://github.com/nextflow-io/nextflow) **3526★**
- 核心逻辑：agentic 的 map-reduce 标注管线，把凌乱的公开**测序样本元数据**（ENA/SRA/GEO）映射为受控词表 ontology 术语。
- **标注对象是样本元数据**（organism / tissue·cell type / disease / assay），**不是序列本身**。
- 记录理由：**是**——`SampleLabels` 含 `rationale`（引用触发判定的原文短语）+ `confidence`(0–1) + 低置信进 `review_queue`。
- **map-reduce 骨架可借鉴**：`META_FETCH(shard) → classify(map) → collect(gather) → audit(reduce)`，把"map 粒度"从"样本"换成"序列位置"，把"skill 词表"换成 IUPAC/修饰语义词表即可迁移到序列层。

### 7. bio-rs（Rust）
- 仓库 [bio-rs/bio-rs](https://github.com/bio-rs/bio-rs) · PyPI [biors](https://pypi.org/project/biors/) · **3★**
- 核心逻辑：AI-ready 的生物数据 I/O、校验、tokenization 引擎。
- **`--profile dna-iupac`**：是**校验 + tokenization 配置档**，不是语义清理。合法模糊 IUPAC 符号→unknown token + `warning(ambiguous_symbol)`；非法符号→unknown token + `error(invalid_symbol)`；保持位置对齐、不重写碱基。
- 逐位判断：**是**（逐位 W/E 标记、结构化 issue code + message + recovery_hint + provenance 哈希）；但这是**诊断理由**，非"这个模糊碱基在生物学上意味着什么"的语义归因；**不做重写/清理决策**。

### 8. atomworks.io（RosettaCommons/atomworks）
- 仓库 [RosettaCommons/atomworks](https://github.com/RosettaCommons/atomworks) · **465★**
- 核心逻辑：`atomworks.io` 是**结构文件**（PDB/mmCIF/SMILES）的化学清理/标准化（去离去基团、修键级/电荷、多构象处理）。
- **对非标准"序列字符"不做清理**——序列层只有 FASTA 的格式解析（括号记法、非天然氨基酸、DNA/RNA 消歧），无语义分类、无理由、无语义感知。

### 9. Data-Juicer `special_characters_filter`（阿里）
- 仓库 [datajuicer/data-juicer](https://github.com/datajuicer/data-juicer) · 源码 [special_characters_filter](https://datajuicer.github.io/data-juicer/en/v1.5.3/_modules/data_juicer/ops/filter/special_characters_filter.html) · **7245★**
- 核心逻辑：`special_char_ratio = len([c for c in text if c in SPECIAL_CHARACTERS]) / len(text)`，比例落在 `[min_ratio=0.0, max_ratio=0.25]` 内才保留，超阈值**整条删除**。
- **基于"比例阈值"，纯字符集合成员计数，无逐位语义判断**。它是通用 LLM 文本清洗工具，根本不理解 IUPAC——且 R/Y/N 等是字母，甚至不被计入"特殊字符"。

### 10. 关键词检索（漏网之鱼排查）
- `semantic sequence cleaning agent` / `IUPAC ambiguity aware cleaning` / `biological data curation agent sequence` / `sequence quality curation LLM agent`：**未找到直接竞品**。
- 最接近的相邻工作均在"变异/阈值/元数据"层，非"序列字符逐位语义"层：BioPerl `Bio::Tools::IUPAC`（模糊码**展开**为枚举）、minSNPs `resolve_IUPAC_missing`（random/most_common **替换**）、DFAST `replace_ambiguous.py`（模糊码全部 **→N**）、BaseCamp / AVA（agentic 决策，但工作在 read 级 variant/VAF 阈值，非字符层）。

## 二、对比表

| 工具 | 对非标准字符的处理策略 | 逐位判断 | 记录理由 | 按下游任务标注 | 语义感知 | ★ |
|---|---|---|---|---|---|---|
| seismic cleanfa | 模糊码→N 替换；非 IUPAC 报错 | 逐位(哑替换) | 否 | 否 | 否 | 14 |
| SNIP | 含非 20AA → 整条删除 | 检测逐位/动作整条 | 部分(计数) | 否 | 否 | 0 |
| ViralEntropR | B/J/X/Z → 整行删除；保留 gap | 检测逐位/动作整行 | 部分(行级) | 否 | 部分(gap 二分类) | 0 |
| ShortRead clean() | 含非 ACGT IUPAC → 整条删除 | 检测逐位/动作整条 | 否 | 否 | 否 | 8 |
| dsh-bio-genie | IUPAC 仅用于序列类型判断 | 否 | 否 | 否 | 弱(字母表归属) | 4 |
| data-labelling | 标注样本**元数据**(非序列) | 逐样本 | **是**(rationale+conf+review) | 是(ontology) | 是(元数据层) | 3526 |
| bio-rs | dna-iupac = 校验+tokenization | **是**(逐位 W/E) | 部分(issue code) | 否 | 部分(符号级) | 3 |
| atomworks.io | 结构层化学清理 | 否 | 否 | 否 | 否 | 465 |
| Data-Juicer | 特殊字符比例阈值 → 整条删 | 否 | 否 | 否 | 否 | 7245 |

## 三、精确空白（回答核心问题）

**问题 A：现有工具把非标准字符当"噪声"还是"信息载体"？**
→ 除"展开/打分"类（BioPerl IUPAC、novoAlign）外，绝大多数把非标准字符**当需要删除/替换/按比例过滤的噪声**：cleanfa/DFAST/minSNPs 替换为 N 或随机碱基，SNIP/ShortRead/ViralEntropR 见一个模糊码就删整条，Data-Juicer 超比例整条丢弃。没有任何工具把它当"需要按生物学上下文分类处理的信息载体"。

**问题 B：它们能区分"测序乱码"和"有生物学意义的 IUPAC 模糊码"吗？**
→ **不能**。现有逻辑是全局/集合级的：要么"非 ATCG 一律删/换"，要么"非标点字母数字按比例过滤"。最接近的是 **bio-rs**（区分合法模糊 `ambiguous_symbol` warning vs 非法 `invalid_symbol` error），但它区分的是"合法 vs 非法**符号**"，不是"这个模糊码**在生物学上意味着什么、该保留还是丢弃**"；且它不重写、不输出"为什么这样处理"的语义理由。

**问题 C：是否已有人做"语义感知序列清理 / 逐位判断非标准字符生物学意义"？**
→ **未找到直接竞品**。最接近的 agentic 工作（BaseCamp、AVA）停留在变异/VAF/阈值决策层；data-labelling 虽具备"rationale + confidence + review_queue"的完整可审计结构，但对象是样本元数据。**"序列逐位语义分类 + 上下文/下游任务感知决策 + 可审计理由"三者合一的切口，目前是空白。**

## 四、本 Skill 的原创性定位

现有工具做的是**字符级规则清洗**（哑替换/整条删除），本 Agent 做的是**语义级判断与记录**：

- 把一个非标准字符归入 `sequencing_error / iupac_ambiguity / functional_modification / unknown_artifact` 四类之一，并给出引用具体依据的理由；
- 同一个字符（如 `N`）在不同上下文（保守基序内 vs 低质量同聚物区）得到**不同语义分类与动作**；
- 决策随**下游任务**（预训练/分类）变化，标签与逐位决策绑定；
- 工具查不到就**降级为人工复核**，绝不猜测。

架构上直接**借鉴 data-labelling 的 map-reduce + rationale + confidence + review_queue**，并把 map 粒度从"样本"下沉到"序列位置"；校验层借鉴 **bio-rs 的逐位符号级 issue code**，但在其之上叠加"生物学语义归因 + 决策动作 + 下游标签"。
