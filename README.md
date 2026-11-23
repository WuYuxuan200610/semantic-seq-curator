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

这里有一条不容易想到、但代价很高的坑：**从已损坏数据派生出的注释列是循环证据**。
本项目的结构注释列曾经"恰好"显示占位区域的宽度等于标记文本宽度，差点让人得出错误结论——
它反映的是损坏后的状态，永远不能用来推断损坏前有多少内容。上面这个方法就是在纠正它。

---

## 二、与邻近工具的差别

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
| 未知约定实证 | 不涉及 | 锚定法脚本 + 合成演示样例 |

补充两点实际影响：

- **数据产出**：占位符按计数还原后，cleanfa 与朴素删除脚本会丢失的那部分长度被补了回来
  （实测 100.2% vs 93.5%，见 §3.5）。
- **注释列寄存器保持**：不做同步时实测约 21.7% 的行会序列/注释错位（见 §3.3）。
- **工具用错域会静默清零数据**：把蛋白白名单工具用在核酸数据上会把序列整条删光。
  这类"结果为空但不报错"的情形，正是上面第 4 条防护要拦的。

---

## 三、实测数据

以下全部数字来自**一份本地原始数据**（未随仓库提供，故只给数值与复现命令）。
每组实验都可用你自己的数据复跑。

### 3.1 清理与独立复核的规模

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

### 3.2 变换账本（机器可读、逐行可查）

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

### 3.3 寄存器保持：不做同步会错位多少？（A/B 实验）

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

### 3.4 占位符 `*[k]*`：k 的含义是被实证的，不是猜的

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

### 3.5 产出可用性对比

```bash
python benchmark_model_ready.py <raw.csv> --seq-col Sequence --chart benchmark/model_ready.png
```

![目标导向对比](benchmark/model_ready.png)

| 方法 | 输出合法率 | 数据产出 | 序列保留 |
|---|---|---|---|
| 本项目规则 | 100.00% | **100.2%** | 全部 |
| seismic cleanfa | 100.00% | 93.5% | 全部 |
| 朴素删除脚本 | 100.00% | 93.5% | 全部 |
| SNIP（蛋白白名单） | — | — | **0** |

- **合法率**：输出字符 ∈ 目标字母表 的比例。<100% 就意味着有字符喂不进模型。
- **数据产出**：输出字符数 / 原始字符数。本项目 >100% 是因为占位符被还原为缺失的碱基，
  这部分长度在另外两种输出里是丢失的。
- SNIP 是蛋白工具（白名单仅 20 种标准氨基酸），用在核酸数据上会**整条删光且不报错**——
  这正是 §1.4 那类防护要拦的情形。

另一张图（[benchmark/comparison.png](benchmark/comparison.png)）是在**合成用例**上
对语义判断正确性的对比，与真实数据无关。

---

## 四、快速开始

```bash
# 1) 清理（自动检测序列列；输出保证 100% 在目标字母表内）
python curate.py <sequences.csv> --out cleaned_output

# 2) 独立复核（不复用清理代码，重新逐字符检查产物）
python validate_model_ready.py cleaned_output/cleaned.csv \
    --align-col rna_seq=structure_curated

# 3) 遇到未知的占位约定，先实证再处理
python verify_placeholder_semantics.py <raw.csv> \
    --seq-col Sequence --ref-col FASTA_sequence --flank 15

# 4) 自检
python test_curate_entry.py    # 入口回归（含"全表无字母表外字符"断言）
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

## 五、架构

```
curate.py                       清理主入口（字母表驱动 + 运行内合法性断言）
validate_model_ready.py         独立复核器（刻意不复用清理实现）
verify_placeholder_semantics.py 未知占位约定的实证方法（含合成演示样例）
benchmark_model_ready.py        合法率 / 数据产出 对比 + 画图
benchmark_alignment.py          寄存器保持 A/B 量化
benchmark.py / draw_chart.py    合成用例上的语义判断对比与绘图
semantic_seq_curator/           语义分类内核
  ├─ iupac.py                   IUPAC 映射表、修饰表、保守基序库
  ├─ downstream.py              下游任务 schema（词表容忍度与标签要求）
  ├─ tools.py                   Tool 抽象 + 工具注册表
  ├─ agent.py                   Reasoner（规则 / LLM 可插拔）+ ReAct 循环
  └─ autonomous.py              由 LLM 面对原始序列自主发现语义 + 规则沉淀
test_curate_entry.py            入口回归测试
test_agent.py                   语义分类组件回归测试
run_demo.py                     合成用例演示
docs/                           竞品调研、差异说明、Skill 定义、System Prompt
examples/                       合成数据（可直接跑通全部脚本）
benchmark/                      对比结果与图
```

**分层设计**：规则引擎是主体，可独立完成全部清理（默认路径、零外部依赖、可离线规模化）。
`Reasoner` 是接口，`RuleReasoner`（确定性、可单测）与 `LLMReasoner`（注入任意大模型）
可互换——需要让模型参与长尾判断时替换推理器即可，清理主链路不变。

---

## 六、边界与已知限制

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

## 七、相关工作与致谢

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
