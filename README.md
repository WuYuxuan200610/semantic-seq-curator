# semantic-seq-curator

把生物序列变成**模型可读、可验证、可审计**的训练输入。

一句话定位：**不只是"把乱码换成 N"，而是给"清洗"这件事补上一份可执行的契约**——
先声明目标词表，再逐字符映射，运行内自检，最后用一份**不复用清理代码**的独立脚本复核；
每个字符的变换都作为数据导出，可逐条复查。

---

## 一、为什么是"契约"而不是"又一个替换脚本"

大多数序列清洗的产物是**一条干净的序列**，你无从知道它被改了什么。本项目的产物是
**序列 + 变换账本 + 合法性证明**三件套。具体做成了五件事：

### 1. 合法性是一个被证明的断言，不是口头承诺

清理时逐行断言，产物写完后**再独立复查一遍**：

- `curate.py` 运行内检查每个输出字符是否 ∈ 目标字母表；
- [`validate_model_ready.py`](validate_model_ready.py) **故意重新实现一遍字母表定义**，
  不 import 清理代码，重新读产物逐字符复核，并统计"出现过的禁用字符集合"。

> 验证器与被验证对象共享实现，验证就没有意义。所以这份验证器是**独立副本**。

只要有一个字符不在词表里，它就喂不进模型——这条检查保证这件事不会静默发生。

### 2. 每个字符的变换都可导出、可追查

产物里每个序列列都带一列 `<列名>_curation_summary`（逐行 JSON）和列级映射直方图：

```json
{"column": "rna_seq", "sequence_type": "rna", "target_alphabet": "rna",
 "original_length": 30, "curated_length": 26,
 "mapping": {"小写→大写": 4, "gap→fallback": 1, "占位符→fallback×k": 1}}
```

"哪一列、改了多少个、改成什么"，全部是机器可读的数据，而不是隐藏在实现里。
清理出问题时可定位到具体行、具体映射项。

### 3. 寄存器保持：序列和注释列一起动，不许错位

这是本项目与"单列文本替换"最实质的区别。

序列和它的注释（二级结构、比对串、位置索引）必须逐位对齐。一旦某个变换**会改变长度**
（比如占位符还原成 `k` 个 `N`），只改序列就会让注释整体错位——而这种错位是静默的，
下游不会报错，只会训出坏模型。

所以 `--align-col` 把注释列纳入同一次变换：同一条规则、同一个位置映射，同步改写。
并且：**无法证明两列原本对齐时不会硬猜**——配对失败直接报错停止；
个别行长度不一致则跳过该列并在该行记录原因。

```bash
python curate.py data.csv --align-col Structure              # 自动配对
python curate.py data.csv --align-col rna_seq=structure      # 显式指定配对
```

### 4. 失败即拒绝，不静默猜

清洗最危险的失败不是报错，而是**悄悄猜错**。因此以下情况一律停止并给出可操作的提示：

| 情况 | 行为 |
|---|---|
| 显式指定的列名不存在 | 报错，列出可用字段，**不写出任何结果** |
| 自动检测不到可信序列列 | 报错并提示用 `--seq-col` 指定，**不写出任何结果** |
| 注释列无法与任何序列列配对 | 报错并提示用 `序列列=注释列` 显式配对 |
| 注释列与序列原始长度不一致 | 跳过该列并记录原因，**不冒险做位置映射** |
| 输出出现字母表外字符 | 复核失败，退出码非 0 |

### 5. 未知的数据约定，用实证而不是猜测定下来

上游工具可能留下 `*[k]*` 这类占位标记，而 `k` 的含义**决定了整个处理方式**。
本仓库把"如何定这件事"做成了一个可复用脚本：
[`verify_placeholder_semantics.py`](verify_placeholder_semantics.py)
——取占位符两侧锚点，在**未被污染的参照列**里量出真实跨度，直接判定 `k` 是计数还是编号。

```bash
# 仓库自带合成样例，可直接跑通
python verify_placeholder_semantics.py examples/placeholder_anchor_demo.csv \
    --seq-col seq_with_placeholder --ref-col reference_sequence --flank 15
# → 真实跨度 == k : 3/3 = 100.00%
# → 真实跨度 == 标记宽度 : 3/3 = 0.00%
```

这里有一个容易踩、代价又高的坑：**从已损坏数据派生出的注释列是循环证据**。
如果注释列本身是由**已经损坏**的序列算出来的（例如结构串），它反映的就是损坏后的状态，
用它推断"损坏前有多少内容"会得到自洽但错误的结论。
判断长度类约定必须找**未被污染的同类参照列**——上面这个方法就是为此设计的。

---

## 二、LLM 到底在哪里起作用（以及哪里不需要它）

这个问题必须先答清楚，否则"skill / agent"就只是包装。

**清洗本身不需要模型。** 把字符映射进已知词表是确定的：规则更快、更稳、更便宜、可复现，
所以 `curate.py` 的默认路径**一行 LLM 都不调**，全量跑是零成本。

**模型不可替代的地方只有一处：面对一份陌生的脏数据，判断那些记号是什么意思。**
脚本只能执行人已经想好的规则；没人告诉它 `*[k]*` 是什么时，它无从下手。

### 内置的"语义专家"模块：`discover_conventions.py`

这就是让本仓库配得上 skill 的模块，也是为什么必须引入 LLM。它给模型划了一条**不可越界的线**：

```
LLM 负责「提出可证伪的假设」「提出检索式」「从检索结果里提取可引用结论」
代码负责「跑实验判定真假」「校验引用是否真实存在」
```

| 阶段 | 谁做 | 内容 |
|---|---|---|
| 1 观察 | 代码（确定性） | 扫描数据，找出可疑记号：`*[k]*`、小写、gap、裸数字、模糊码……给出分布和样例 |
| 2 提议 | **LLM** | 对每个记号提出假设，并**只能从受支持的预测清单里选一个**来验证 |
| 3 验证 | 代码（确定性） | 用参照列做侧翼锚定，量出真实跨度；模型不参与判定 |
| 4 外部证据 | **LLM 提检索式并提取结论；代码校验引用** | 数据内判不了时，去查工具文档/论文/流程/论坛——**去查，而不是猜** |
| 5 裁定 | 代码（确定性） | 按证据等级给出最终结论（见下表） |

**关键设计：模型不能说"看起来像"。** 它必须说"我预测 X 成立"，
然后由代码去量。它挑不出验证方案，该假设直接标 `UNTESTABLE`，绝不当成通过。
这既用上了模型对陌生格式的联想能力，又不给它幻觉的机会。

### 证据是分层的：`DOCUMENTED` ≠ `CONFIRMED`

像人类做研究一样，不同来源的证据强度不同，报告必须标明用的是哪一层：

| 等级 | 证据来源 | 强度 |
|---|---|---|
| `CONFIRMED` | **数据内实验**（参照列锚定，可复现） | 最强 |
| `DOCUMENTED` | **外部来源明确记载**：工具官方文档 / 论文 / 分析流程（须给出 URL 与原文片段） | 强，但未经数据验证 |
| `REFUTED` | 数据内实验证伪 | — |
| `WEAK` | 命中率介于 5%~95%，证据不足，**不得视为成立** | 弱 |
| `UNVERIFIED_HYPOTHESIS` | 模型凭记忆提出、无可引用来源 | 仅作线索 |
| `UNTESTABLE` | 无参照列、检索也无命中 | 不猜 |

**防幻觉的硬约束**：模型引用的 URL **必须真的出现在检索结果里**，支撑片段必须与来源内容相符，
否则结论被拒、降级为 `UNTESTABLE`。`test_discover.py` 中有对应断言。

实测（本地原始数据，`--flank 15`，数据内有参照列）：

| 记号 | 假设（LLM 提出） | 预测 | 样本 | 命中率 | 裁定 |
|---|---|---|---|---|---|
| `*[k]*` | 计数 | 跨度 == 记号里的数字 | 3,000 | **100.0%** | **CONFIRMED** |
| `*[k]*` | 固定宽度注解 | 跨度 == 记号文本宽度 | 3,000 | 2.1% | **REFUTED** |
| 小写字母 | 低置信度位点 | 跨度 == 记号宽度（1） | 3,013 | **100.0%** | **CONFIRMED** |

第三条的结论有实际含义：**小写字母恰好占据一个真实碱基位置**，说明它是"低置信度的碱基调用"，
而不是被插入的注解——这正是"转大写并保留"这个处理动作的依据，而它是**被实验证明的，不是被规定的**。

当**没有参照列**时，流程自动转入第 4 阶段去查资料。离线合成演示（mock 检索，
URL 用 `example.invalid` 占位，**不是真实文献**）：

| 记号 | 数据内 | 外部证据 | 最终裁定 |
|---|---|---|---|
| `*[k]*` | UNTESTABLE | 命中官方文档，有原文片段 | `DOCUMENTED` |
| 小写字母 | UNTESTABLE | 未查到关于它的记载 | `UNTESTABLE`（不硬套） |

确认的假设会输出成规则建议，人工过目后才交给 `curate.py`，并标注证据等级：

```json
[{"regex": "\\*\\[ *(\\d+)\\]\\*", "action": "expand_to_n_by_number",
  "evidence_tier": "CONFIRMED", "evidence": "数据内锚定 3000 例，命中率 1.0",
  "requires_small_sample_check": false, "maps_to": "占位符 -> k 个 N"}]
```

### 外部证据：脚本不碰网络，检索交给 agent

**为什么不让脚本自己检索？** 实测过：把搜索结果页抓下来再解析 HTML 非常脆弱——
查询词稍偏就返回完全不相关的内容（搜序列记号却返回视频网站），
中文搜索引擎对英文技术术语召回也很差。真正管用的是**读页面内容、判断相关性、换词重试**，
这是 agent 的强项，不是脚本的强项。

所以这里做成**两阶段的文件交换**，脚本保持纯确定性、不发任何网络请求：

```bash
# 第一趟：脚本产出「待查清单」——要查什么、什么算合格证据、结果怎么写回来
python discover_conventions.py <raw.csv> --seq-col Sequence --emit-queries queries.json

# 中间：由具备联网能力的 agent（如 DSH 会话）读清单、执行检索、写出证据文件
#       agent 用它自己的 web_search / web_fetch，无需本仓库配置任何 API key

# 第二趟：脚本读入证据，做与自动检索完全相同的防幻觉校验，然后裁定
python discover_conventions.py <raw.csv> --seq-col Sequence \
    --evidence evidence.json --emit-rules rules.json
```

清单长这样（节选）：

```json
{"instructions": ["对 unresolved 里的每一项，使用你的联网检索能力查明该记号的约定", "..."],
 "evidence_quality_rules": {"DOCUMENTED": "需要 official_docs / peer_reviewed / protocol 级来源，且 url 与原文片段都要能对应上"},
 "unresolved": [{"pattern": "star_bracket_number", "description": "星号方括号数字记号，如 *[19]*",
                 "occurrences": 6039, "question": "什么工具会产生这种记号？",
                 "suggested_queries": ["... sequence notation meaning", "..."]}],
 "response_schema": {"<pattern>": {"search_queries": [...], "hits": [{...}], "claim": {...}}}}
```

**校验不因来源不同而放松**：即使证据是 agent 提供的，URL 仍必须真实出现在该记号的 hits 里、
支撑片段仍必须与来源相符，否则结论被拒、降级 `UNTESTABLE`。

**检索顺序由待查清单给出**：代码搜索优先——这类记号的权威定义通常就写在生成它的源码里；
其次是产生该记号的工具的官方文档，再次论文/分析流程文档，最后论坛。

agent 用它的任一联网工具完成检索，把结果写成证据文件即可，**脚本无需改动**。

也保留一条**全自动路径**（`--search-cmd`），适合已有搜索引擎 API 的场景：

```bash
# 检索桥：从 stdin 读查询串 → 把 JSON 结果列表写 stdout，每条 {title, url, snippet}
python discover_conventions.py <raw.csv> --seq-col Sequence \
    --llm-cmd "python llm.py" --search-cmd "python search.py" --emit-rules rules.json
```

> ⚠️ **验证状态**：两阶段文件交换与防幻觉校验已完整跑通并有单测覆盖
> （含 BOM 容忍、字段名兼容、编造引用被拒），离线 mock 检索也已跑通。
> 但**接上真实搜索引擎后的检索质量未做端到端验证**。

**离线可复现**：仓库自带 `examples/convention_discovery_demo.csv`（合成，60 行）：

```bash
# 有参照列 → 走数据内验证
python discover_conventions.py examples/convention_discovery_demo.csv \
    --seq-col seq_with_markers --ref-col reference_sequence --llm mock

# 无参照列 → 产出待查清单（离线可跑；真实检索需 agent 参与）
python discover_conventions.py examples/convention_discovery_demo.csv \
    --seq-col seq_with_markers --llm mock --emit-queries queries.json

# 无参照列 + 离线合成检索（仅用于跑通全流程）
python discover_conventions.py examples/convention_discovery_demo.csv \
    --seq-col seq_with_markers --llm mock --search mock
```

**换一份陌生数据，这个循环照样跑**：观察 → 提假设 → 锚定验证 → 查资料 → 裁定。
这就是"下载后可泛化使用"的那部分。

---

## 三、与邻近工具的差别

这些工具都很成熟，本项目不是它们的替代品，而是在它们没有覆盖的位置上工作：

| 能力 | seqkit / cleanfa / SNIP 等 | 本项目 |
|---|---|---|
| 字母表规范化 | ✅ 成熟、快、生态好 | ✅ 同样支持 |
| 目标词表可切换 | cleanfa 面向核酸、SNIP 面向蛋白，各自固定 | **同一入口覆盖核酸 + 蛋白**，`dna` / `rna` / `nucleic-mixed` / `dna-iupac` / `protein` / `protein-extended` 随 tokenizer 切换 |
| 词表支持模糊码时保留信息 | 通常不具备 | `--alphabet dna-iupac` **原样保留** R/Y/W/… 而不是丢弃 |
| 变换账本（改了什么、为什么） | 一般只给结果或计数 | **逐行 JSON + 列级直方图**，可导出复查 |
| 独立复核脚本 | 罕见 | **刻意不复用实现**的验证器 |
| 注释列寄存器保持 | 不在其职责内 | `--align-col` 同步变换；无法证明对齐则拒绝 |
| 静默猜错防护 | 视工具而定 | 列/类型/配对不确定一律停止，**不写结果** |
| 未知约定实证 | 不涉及 | **`discover_conventions.py`：LLM 提假设、代码跑实验判定**（见 §二），可在任何陌生数据上泛化 |
| 可泛化的语义专家模块 | 无（规则写死） | 观察 → 提假设 → 锚定验证 → 裁定的完整循环，含受支持的预测清单与 `UNTESTABLE` 保守兜底 |

补充两点实际影响：

- **数据产出**：占位符按计数还原后，cleanfa 与朴素删除脚本会丢失的那部分长度被补了回来
  （实测字符净变化 **+8,731** vs **−367,530**，见 §4.5）。
- **注释列寄存器保持**：不做同步时实测约 21.7% 的行会序列/注释错位（见 §4.3）。
- **工具用错域会静默清零数据**：把蛋白白名单工具用在核酸数据上会把序列整条删光。
  这类"结果为空但不报错"的情形，正是上面第 4 条防护要拦的。

---

## 四、实测数据

以下全部数字来自**一份本地原始数据**（未随仓库提供，故只给数值与复现命令）。
每组实验都可用你自己的数据复跑。

### 4.1 清理与独立复核的规模

`curate.py` 自检 + `validate_model_ready.py`（不复用实现的独立副本）复核：

| 指标 | 数值 |
|---|---|
| 行数 / 序列列数 | 23,535 / 8 |
| 复核字符数 | **23,509,650** |
| 字母表外字符（独立复核） | **0** |
| 出现过的禁用字符集合 | `[]`（空） |
| 输出实际出现的字符全集 | `ACGNTU` |
| 序列与注释列错位的行 | **0** |

```bash
python curate.py <raw.csv> --align-col Structure --out out
python validate_model_ready.py out/cleaned.csv --align-col Sequence=Structure_curated
```

### 4.2 变换账本（机器可读、逐行可查）

同一份数据上，每一个被改动的字符都被计数，并可下钻到具体行：

| 映射项 | 数量 |
|---|---|
| 小写 → 大写 | 524,997 |
| gap → `N` | 327,843 |
| 占位符 `*[k]*` → `k` 个 `N` | 6,680 个标记 → **48,418** 个字符 |
| IUPAC 模糊码 → `N` | 2,097 |
| 长度发生变化的行 | 5,109 |

这正是"清洗"变成"契约"的地方：**改了多少、改成什么，是产物的一部分**，
而不是隐藏在实现里。逐行明细在每行的 `<列名>_curation_summary`，列级汇总在同名 JSON。

### 4.3 寄存器保持：不做同步会错位多少？（A/B 实验）

同一份数据跑两次，唯一差别是加不加 `--align-col`：

| 配置 | 行数 | 序列/注释错位行 | 错位率 |
|---|---|---|---|
| **A** `--align-col`（注释列同步变换） | 23,534 | **0** | **0.000%** |
| **B** 不加（注释列保持原样） | 23,534 | 5,109 | **21.709%** |

**不做同步，约 21.7% 的行会出现序列与注释错位**——而且这种错位不会报错，
下游只会训出坏模型。这是本项目最实质的一项收益，也是它区别于"单列文本替换"的地方。

```bash
python benchmark_alignment.py <raw.csv> --align-col Structure
# 仓库自带合成样例可离线复现（3 行，B 配置 100% 错位、A 配置 0）
python benchmark_alignment.py examples/example_sequences.csv \
    --seq-col rna_seq --align-col structure
```

### 4.4 占位符 `*[k]*`：k 的含义是被实证的，不是猜的

用未被污染的参照列做侧翼锚定（`--flank 15`），量出真实跨度：

| 判据 | 命中 / 可测样本 | 比例 |
|---|---|---|
| 真实跨度 == `k`（计数模型） | 2,168 / 2,168 | **100.00%** |
| 真实跨度 == 标记文本宽度（编号模型） | 36 / 2,168 | 1.66% |

结论：`k` 是该处丢失的碱基个数 → 应还原为 `k` 个 `N`。

> 那 1.66% 的残余不是矛盾：当 `k` 与标记宽度**恰好相等**时（如 `*[5]*` 宽 5），
> 两个判据本来就会给出同一答案。有判别力的样本上，两个模型完全分离。
>
> 另一个必须记住的教训：锚点太短会误匹配。`--flank 10` 时出现过 3/2967 的假阳性，
> 加长到 15 后消失。**从已损坏数据派生出的注释列（如结构串）是循环证据**，
> 不能用来自证长度——必须找未被污染的参照。

```bash
python verify_placeholder_semantics.py <raw.csv> \
    --seq-col Sequence --ref-col FASTA_sequence --flank 15
# 离线自检（合成样例）
python verify_placeholder_semantics.py examples/placeholder_anchor_demo.csv \
    --seq-col seq_with_placeholder --ref-col reference_sequence --flank 15
```

### 4.5 产出可用性对比

```bash
python benchmark_model_ready.py <raw.csv> --seq-col Sequence --chart benchmark/model_ready.png
```

![目标导向对比](benchmark/model_ready.png)

| 方法 | 输出合法率 | 输出字符数 | 字符净变化 | 序列保留 |
|---|---|---|---|---|
| 本项目规则 | 100.00% | 5,652,399 | **+8,731** | 全部 |
| seismic cleanfa | 100.00% | 5,276,138 | −367,530 | 全部 |
| 朴素删除脚本 | 100.00% | 5,275,980 | −367,688 | 全部 |
| SNIP（蛋白白名单） | — | 0 | −5,643,668 | **0** |

**图里只画了前三家**：SNIP 的输出是"整条清零"，属性质不同的结果，不是同一量纲上的极端值；
把它画进同一坐标会把另外三根柱子压成不可见。它的结果由上表与图下脚注说明。

- **合法率**：输出字符 ∈ 目标字母表 的比例。<100% 就意味着有字符喂不进模型。
  输出 0 个字符时合法率无从计算，记作「—」。
- **净变化**：输出字符数 − 输入字符数。**这里用绝对字符数而不是百分比**，
  因为比值会超过 100%，容易被误读成"注水"。
  `>0` 的含义是：**输入本身已经被上游删掉了一部分碱基**（占位符就是那些缺失的痕迹），
  还原后比输入长。而 cleanfa / 朴素脚本会把这段长度连同 gap 一起丢掉，所以是负数。

  换句话说，本项目在这份数据上比"输入"多补回 8,731 个字符，
  而其他两种做法比输入少 36.7 万个——差距不只在百分比，而在**能否对齐原始生物学长度**。
- SNIP 是蛋白工具（白名单仅 20 种标准氨基酸），用在核酸数据上会**整条删光且不报错**——
  这正是 §1.4 那类防护要拦的情形。

另一张图（[benchmark/comparison.png](benchmark/comparison.png)）是在**合成用例**上
对语义判断正确性的对比，与真实数据无关。

---

## 五、快速开始

### 获取

```bash
git clone https://github.com/WuYuxuan200610/semantic-seq-curator.git
cd semantic-seq-curator
```

只需要 Python 3.9+，**无第三方依赖**（`matplotlib` 仅画图时需要，可选）。

### 跑一遍

```bash
# 1) 清理（自动检测序列列；输出保证 100% 在目标字母表内）
python curate.py <sequences.csv> --out cleaned_output

# 2) 独立复核（不复用清理代码，重新逐字符检查产物）
python validate_model_ready.py cleaned_output/cleaned.csv \
    --align-col rna_seq=structure_curated

# 3) 遇到陌生记号：让模型提假设、由代码跑实验判定（离线可跑）
python discover_conventions.py examples/convention_discovery_demo.csv \
    --seq-col seq_with_markers --ref-col reference_sequence --llm mock

# 3b) 没有参照列时：产出待查清单 → agent 联网检索 → 喂回证据裁定
python discover_conventions.py <raw.csv> --seq-col Sequence --emit-queries queries.json
#    （agent 按清单检索后写出 evidence.json）
python discover_conventions.py <raw.csv> --seq-col Sequence --evidence evidence.json \
    --emit-rules rules.json

# 4) 寄存器保持量化（同一数据开/关 --align-col 跑两次）
python benchmark_alignment.py <raw.csv> --align-col Structure

# 5) 自检
python test_curate_entry.py    # 入口回归（含"全表无字母表外字符"断言）
python test_discover.py        # 发现循环回归（观察→假设→验证→裁定）
python test_agent.py           # 语义分类组件回归
python run_demo.py             # 合成用例：同一字符在不同上下文得到不同决策
```

依赖：仅标准库。画图需要 `matplotlib`（可选）。

### 目标字母表

| 档案 | 允许字符 | 适用 |
|---|---|---|
| `dna` | A C G T N | DNA 语言模型 |
| `rna` | A C G U N | RNA 语言模型 |
| `nucleic-mixed` | A C G T U N | 不能确定 DNA/RNA 时 |
| `dna-iupac` / `rna-iupac` | 标准碱基 + 全部 IUPAC 模糊码 + N | tokenizer 支持模糊码时**保留信息** |
| `protein` | 20 种标准氨基酸 + X | 严格蛋白词表 |
| `protein-extended` | 20 + X + B Z J U O * | ESM 类词表 |
| `auto`（默认） | 按行推断类型自动选择 | 通用 |

### 映射规则

| 输入 | 处理 |
|---|---|
| 小写 `acgu` | → 大写 |
| IUPAC 模糊码 `R Y W S K M …` | → `N`；`-iupac` 档案则原样保留 |
| gap `-` `.` | → `N`（默认，保持对齐）或 `--gaps drop` |
| 占位符 `*[k]*` | → `k` 个 `N`（前提：已用锚定法确认 `k` 的含义） |
| 其它未知字符 | → `N`/`X`（默认）或 `--unknown drop` |

### 作为库使用

```python
from curate import clean_sequence
cleaned, mapping, seq_type, alphabet = clean_sequence("ACGTNacgt-*[ 2]*", alphabet="dna")
# cleaned = "ACGTNACGTNNN"     ← 每个字符都是合法 token
# mapping = {"小写→大写": 4, "gap→fallback": 1, "占位符→fallback×k": 1}
```

---

## 六、架构

```
curate.py                       清理主入口（字母表驱动 + 运行内合法性断言）
discover_conventions.py         语义专家模块：LLM 提假设 + 代码跑实验判定（§二）
validate_model_ready.py         独立复核器（刻意不复用清理实现）
verify_placeholder_semantics.py 单一实验的独立实现（锚定法，便于单独复跑）
benchmark_model_ready.py        合法率 / 字符净变化 对比 + 画图
benchmark_alignment.py          寄存器保持 A/B 量化
benchmark.py / draw_chart.py    合成用例上的语义判断对比与绘图
semantic_seq_curator/           逐位语义分类内核
  ├─ iupac.py                   IUPAC 映射表、修饰表、保守基序库
  ├─ downstream.py              下游任务 schema（词表容忍度与标签要求）
  ├─ tools.py                   Tool 抽象 + 工具注册表
  ├─ agent.py                   Reasoner（规则 / LLM 可插拔）+ ReAct 循环
  └─ autonomous.py              让 LLM 直接对序列做语义分类的接口
test_curate_entry.py            入口回归测试
test_discover.py                发现循环回归测试
test_agent.py                   语义分类组件回归测试
run_demo.py                     合成用例演示
docs/                           竞品调研、差异说明、Skill 定义、System Prompt
examples/                       合成数据（可直接跑通全部脚本）
benchmark/                      对比结果与图
```

**两层分工**（这也是"为什么是 skill"的答案）：

| 层 | 文件 | 是否需要 LLM | 职责 |
|---|---|---|---|
| 执行层 | `curate.py` | **否** | 已确认的规则，全量、离线、零成本、可复现 |
| 发现层 | `discover_conventions.py` | **是** | 面对陌生数据提出并验证假设，产出规则 |
| 复核层 | `validate_model_ready.py` | 否 | 不复用实现的独立检查 |

规则引擎（`semantic_seq_curator/`）是执行层的内核：默认路径 `RuleReasoner` 确定性、可单测；
`LLMReasoner` 是可替换接口，需要模型参与长尾判断时替换推理器即可，主链路不变。

---

## 七、边界与已知限制

写清楚边界，比假装没有边界更有用：

1. **约定必须按数据集确认。** 占位符 / gap / 大小写的含义是**逐数据集**的约定，
   换数据请先跑 `verify_placeholder_semantics.py`，不要把某个数据集的结论当通用规律。
2. **占位符只能还原长度，不能还原身份。** `*[k]*` → `k` 个 `N`，原来的碱基是什么**不可知**；
   `N` 表达"未知"，不是"已知是 N"。
3. **`gap → N` 是编码选择，不是生物学结论。** gap 表示相对比对缺失。若要真正的去比对形式，
   用 `--gaps drop`（长度会变，注释列需重新对齐）。
4. **尚未做跨数据集泛化验证。** 现有验证集中于合成用例与一份本地数据。
5. **`LLMReasoner` / `autonomous.py` 需要外部模型后端**才能发挥；默认路径不依赖它们。

---

## 八、相关工作与致谢

| 工具 | 做什么 |
|---|---|
| [seqkit](https://github.com/shenwei356/seqkit) | 序列工具集，覆盖格式转换与字母表规范化 |
| [seismic cleanfa](https://github.com/rouskinlab/seismic-rna) | FASTA 规范：模糊码 → N |
| [Biopython](https://github.com/biopython/biopython) | 序列读写与转换基础库 |
| [ShortRead](https://github.com/Bioconductor/ShortRead) | 短读过滤（Bioconductor） |
| [SNIP](https://github.com/BBCMdP/SNIP) | 按 20 种标准氨基酸白名单整条剔除 |
| [bio-rs](https://github.com/bio-rs/bio-rs) | 逐位符号级校验与 tokenization |
| [data-labelling](https://github.com/nextflow-io/nextflow/tree/master/examples/agents/14_data-labelling) | agentic 元数据标注（rationale + confidence + 复核队列） |

借鉴来源：`data-labelling` 的可审计标注结构；`bio-rs` 的逐位符号级诊断思路；
`dsh-bio-genie` 的序列类型自动判断思路；IUPAC-IUB 标准（Cornish-Bowden 1985, *NAR* 13:3021）。
详细的竞品行为核实见 [docs/RESEARCH.md](docs/RESEARCH.md)。

## License

[MIT](LICENSE)
