"""未知数据约定的发现模块 —— 这是本仓库唯一真正需要 LLM 的地方。

为什么需要 LLM：
    清洗本身（把字符映射进词表）不需要模型，规则更快更稳更便宜——这一点本仓库从不含糊。
    真正需要模型的是**面对一份陌生的脏数据，判断它的记号是什么意思**。
    脚本只能执行人已经想好的规则；没人告诉它 `*[k]*` 是什么时，它无能为力。

但也不能让模型说了算（会幻觉）。所以这里定下一个不可让步的分工：

    LLM 负责"提出可证伪的假设"  →  代码负责"跑实验判定真假"

    ┌─ 阶段 1  观察（确定性）：扫描数据，找出可疑记号及其分布
    ├─ 阶段 2  提议（LLM）  ：对每个记号提出假设，并从**受支持预测清单**中选一个
    ├─ 阶段 3  验证（确定性）：用参照列做锚定等实验，量出真实结果
    ├─ 阶段 4  外部证据     ：数据内判不了时，去查工具文档/论文/流程/论坛（可引用来源），
    │                       而不是猜。模型只能引用检索结果里真实存在的链接。
    └─ 阶段 5  裁定（确定性）：CONFIRMED / DOCUMENTED / REFUTED / WEAK /
                              UNVERIFIED_HYPOTHESIS / UNTESTABLE

LLM 说不出"用什么实验验证"时，该假设直接标为 UNTESTABLE —— 不允许"看起来像"就下结论。

用法：
    # 离线自检（内置 mock 假设，不连任何模型）
    python discover_conventions.py <raw.csv> --seq-col Sequence --ref-col FASTA_sequence --llm mock

    # 接入任意模型：命令从 stdin 读 prompt、把回复写到 stdout
    python discover_conventions.py <raw.csv> --seq-col Sequence --ref-col FASTA_sequence \
        --llm-cmd "python my_llm_bridge.py"

    # 把确认的结论输出成规则建议，人工确认后再给 curate.py 用
    python discover_conventions.py <raw.csv> ... --emit-rules rules.json
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import os
import re
import subprocess
import sys

# ---------------------------------------------------------------------------
# 支持的可证伪预测：LLM 只能从这个清单里选，避免它编造无法执行的验证方案
# ---------------------------------------------------------------------------
SUPPORTED_PREDICTIONS = {
    "span_equals_token_width":
        "锚定出的真实跨度 == 该记号的文本宽度（说明记号只是注解，不占真实碱基）",
    "span_equals_number_in_token":
        "锚定出的真实跨度 == 记号里的数字（说明数字是被替换掉的碱基数）",
    "span_equals_one":
        "锚定出的真实跨度 == 1（说明该处恰好对应一个真实碱基）",
    "span_is_zero":
        "锚定出的真实跨度 == 0（说明该处没有真实碱基）",
    "span_exceeds_token_width":
        "真实跨度 > 记号文本宽度（说明记号展开后应比自身更长）",
}

CONFIRM_THRESHOLD = 0.95
REFUTE_THRESHOLD = 0.05
MIN_SAMPLES = 20          # 少于这个数不足以判定，标 UNTESTABLE

# ---------------------------------------------------------------------------
# 证据分层：像人类生物学家一样，证据有强弱，且必须标明来源
#
#   第 1 层  数据内验证（锚定实验）      → CONFIRMED / REFUTED   —— 最强，可复现
#   第 2 层  外部文献/文档/论坛（可引用） → DOCUMENTED            —— 强，但依赖来源可靠性
#   第 3 层  模型自身记忆（无来源）      → UNVERIFIED_HYPOTHESIS  —— 最弱，只作线索
#
# 关键："去查"不是"猜"。查到的东西带 URL 与出处，可被人复核；
# 模型凭记忆说的一句话没有来源，只能算线索，绝不能当结论。
# ---------------------------------------------------------------------------
SOURCE_TIERS = {
    "official_docs": 1,   # 产生该记号的那个工具自己的文档
    "peer_reviewed": 2,   # 论文
    "protocol": 3,        # 实验/分析流程文档
    "forum": 4,           # 论坛、issue、问答
    "unknown": 9,
}
# 只有这些来源类型足以支撑 DOCUMENTED
DOCUMENTED_SOURCES = {"official_docs", "peer_reviewed", "protocol"}

VERDICT_MEANING = {
    "CONFIRMED": "数据内实验验证成立（最强，可复现）",
    "DOCUMENTED": "外部来源明确记载（有出处可复核，但未经数据验证）",
    "REFUTED": "数据内实验证伪",
    "WEAK": "证据不足，不得视为成立",
    "UNVERIFIED_HYPOTHESIS": "模型凭记忆提出、无可引用来源，仅作线索",
    "UNTESTABLE": "无参照列且未找到外部来源，无法判定（不猜）",
}


# ---------------------------------------------------------------------------
# 阶段 1：观察（确定性，不需要模型）
# ---------------------------------------------------------------------------
SUSPECT_PATTERNS = [
    ("star_bracket_number", re.compile(r"\*\[ *(\d+)\]\*"),
     "星号方括号数字记号，如 *[19]*"),
    ("lowercase_letters", re.compile(r"[a-z]"),
     "小写字母（比对惯例中常代表低置信度/可变位点）"),
    ("gap_characters", re.compile(r"[-.]"),
     "比对缺口字符 - 或 ."),
    ("digit_run", re.compile(r"\d+"),
     "裸露的数字串"),
    ("ambiguous_code", re.compile(r"[NRYWSKMBDHVX]"),
     "IUPAC 模糊码"),
    ("non_alnum", re.compile(r"[^A-Za-z0-9]"),
     "其它非字母数字符号"),
]


def observe(csv_path, seq_col, sample_limit):
    """扫描数据，返回 {列: {字符直方图, 命中的可疑记号}}。"""
    profile = {}
    with open(csv_path, encoding="utf-8-sig", errors="replace", newline="") as f:
        reader = csv.DictReader(f)
        if seq_col not in (reader.fieldnames or []):
            raise SystemExit(f"[错误] 列不存在：{seq_col}")
        hist = collections.Counter()
        found = collections.OrderedDict()
        examples = {}
        rows = 0
        for row in reader:
            value = row.get(seq_col) or ""
            if not value:
                continue
            rows += 1
            hist.update(value)
            for name, regex, desc in SUSPECT_PATTERNS:
                matches = regex.findall(value) if regex.groups else regex.findall(value)
                n = len(matches)
                if not n:
                    continue
                if name not in found:
                    found[name] = {"pattern": name, "description": desc,
                                   "regex": regex.pattern, "count": 0, "distinct": set()}
                found[name]["count"] += n
                for m in matches[:5]:
                    token = m if isinstance(m, str) else m[0]
                    found[name]["distinct"].add(token[:12])
                if name not in examples:
                    hit = regex.search(value)
                    lo = max(0, hit.start() - 25)
                    examples[name] = value[lo:hit.end() + 25]
            if rows >= sample_limit:
                break

    suspects = []
    for name, info in found.items():
        info["distinct"] = sorted(info["distinct"])[:8]
        info["example"] = examples.get(name, "")
        suspects.append(info)
    profile[seq_col] = {
        "rows_scanned": rows,
        "alphabet": "".join(sorted(hist)),
        "char_histogram": dict(hist.most_common(15)),
        "suspects": suspects,
    }
    return profile


# ---------------------------------------------------------------------------
# 阶段 2：提议（LLM）
# ---------------------------------------------------------------------------
PROMPT_TEMPLATE = """你是资深生物信息数据管家。下面是一份序列数据里自动扫描出的可疑记号。
请对每个记号提出一个关于"它到底是什么"的假设，并指定用哪个预测来验证它。

只能使用下列预测之一（必须原样照抄字串）：
{predictions}

可用的验证手段：
- 若提供了参照列（未被污染的同类序列），将用侧翼锚定量出记号的真实跨度，再核对上面选定的预测。

要求：
1. 每个记号给一条假设；不确定就给多条不同假设（会被分别验证）。
2. 不许编造预测名；不许给无法验证的假设（那种情况请把 prediction 留空字符串 ""）。
3. rationale 说明你为什么这么猜，一句话。

只输出 JSON 数组，元素格式：
{{"pattern": "<记号名>", "hypothesis": "<计数|固定宽度注解|低置信度位点|修饰|噪声|其它>",
  "rationale": "<一句话>", "prediction": "<上面清单里的字串，或空串>"}}

扫描结果：
{profile}
"""


def _mock_llm(prompt: str) -> str:
    """离线自检用：按提示词类型分派，模拟真实模型在各阶段的输出。

    真实模式下每次调用是独立的模型请求；mock 需要自己辨认"现在被要求做什么"。
    """
    # 阶段 4a：要求给出检索式
    if "用于检索的查询串" in prompt:
        return json.dumps({"queries": [
            "序列工具 括号数字记号 *[N]* 含义",
            "什么流程会把残基折叠成 *[k]* 记号",
            "bracket number marker sequence alignment convention",
        ]}, ensure_ascii=False)
    # 阶段 4c：要求从检索结果里提取可引用结论
    if "有没有**明确记载**" in prompt or "supporting_snippet" in prompt:
        return _mock_claim(prompt)
    # 阶段 2：要求给出假设
    return json.dumps([
        {"pattern": "star_bracket_number", "hypothesis": "计数",
         "rationale": "数字很像被替换掉的碱基个数", "prediction": "span_equals_number_in_token"},
        {"pattern": "star_bracket_number", "hypothesis": "固定宽度注解",
         "rationale": "它也可能只是个不占真实位置的标记", "prediction": "span_equals_token_width"},
        {"pattern": "lowercase_letters", "hypothesis": "低置信度位点",
         "rationale": "比对惯例中小写是低置信度", "prediction": "span_equals_token_width"},
    ], ensure_ascii=False)


def call_llm(mode, cmd, prompt):
    if mode == "mock":
        return _mock_llm(prompt)
    if not cmd:
        raise SystemExit("[错误] 需要 --llm mock 或 --llm-cmd '<命令>'")
    try:
        result = subprocess.run(cmd, shell=True, input=prompt, capture_output=True,
                                text=True, encoding="utf-8", errors="replace")
    except OSError as exc:
        raise SystemExit(f"[错误] 调用模型失败：{exc}")
    if result.returncode != 0:
        raise SystemExit(f"[错误] 模型命令返回 {result.returncode}：{result.stderr[:400]}")
    return result.stdout


def parse_hypotheses(raw):
    text = raw.strip()
    if text.startswith("```"):
        text = re.sub(r"^```[a-zA-Z]*\n?|```$", "", text).strip()
    try:
        data = json.loads(text)
    except json.JSONDecodeError:
        return [], "模型输出不是合法 JSON"
    if not isinstance(data, list):
        return [], "模型输出不是数组"
    cleaned = []
    for item in data:
        if not isinstance(item, dict):
            continue
        pred = item.get("prediction") or ""
        if pred and pred not in SUPPORTED_PREDICTIONS:
            item["prediction"] = ""
            item["rationale"] = (item.get("rationale", "") +
                                 " [原预测名不受支持，已置空]").strip()
        cleaned.append(item)
    return cleaned, ""


# ---------------------------------------------------------------------------
# 阶段 3：验证（确定性探测器）
# ---------------------------------------------------------------------------
def rna_to_dna(text):
    return text.replace("U", "T").replace("u", "t")


LETTER_RE = re.compile(r"[A-Za-z]")
MIN_FLANK = 8   # 侧翼至少要有这么多字母才够定位


def _flank(seq, start, end, direction, want):
    """取目标记号旁的侧翼，**只收字母**并在收够 want 个后停止。

    只收字母是刻意的：参照列是干净的序列，只有字母能在里面定位；
    一旦把 `*`、`[`、`-` 这些噪声收进侧翼，锚点就永远找不到——
    这会静默地把样本数变成 0（本模块早期的真实缺陷）。
    """
    out = []
    i = start - 1 if direction < 0 else end
    while 0 <= i < len(seq) and len(out) < want:
        if LETTER_RE.match(seq[i]):
            out.append(seq[i])
            i += direction
        else:
            break          # 碰到噪声就停，宁可侧翼短也不污染锚点
    if len(out) < MIN_FLANK:
        return None
    if direction < 0:
        out.reverse()      # 向左收集是逆序的，必须翻回来才能当锚点
    return rna_to_dna("".join(out)).upper()


def probe_anchor_span(csv_path, seq_col, ref_col, regex_src, flank, sample_limit):
    """用参照列做侧翼锚定，量出每个记号的真实跨度。

    返回 {"tested": n, "spans": [...], "records": [(token, span), ...]}
    records 同时用于需要 token 的预测（宽度、数字），因此一次扫描即可。
    """
    regex = re.compile(regex_src)
    spans = []
    records = []
    with open(csv_path, encoding="utf-8-sig", errors="replace", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            seq = row.get(seq_col) or ""
            ref = (row.get(ref_col) or "").upper()
            if not seq or not ref:
                continue
            for m in regex.finditer(seq):
                left = _flank(seq, m.start(), m.end(), -1, flank)
                right = _flank(seq, m.start(), m.end(), +1, flank)
                if left is None or right is None:
                    continue
                li = ref.find(left)
                if li < 0:
                    continue
                ri = ref.find(right, li + len(left))
                if ri < 0:
                    continue
                span = ri - (li + len(left))
                spans.append(span)
                records.append((m.group(0), span))
            if len(spans) >= sample_limit:
                break
    return {"tested": len(spans), "spans": spans, "records": records}


def evaluate_predictions(probe_result):
    """对每个受支持的预测算出命中率。返回 {预测名: 命中率 或 None}。"""
    spans = probe_result["spans"]
    records = probe_result["records"]
    if not spans:
        return {}
    n = len(spans)
    out = {}

    hit_w = sum(1 for token, span in records if span == len(token))
    out["span_equals_token_width"] = hit_w / len(records)

    numeric = [(token, span) for token, span in records if re.findall(r"\d+", token)]
    if numeric:
        hit_num = sum(1 for token, span in numeric
                      if span == int(re.findall(r"\d+", token)[0]))
        out["span_equals_number_in_token"] = hit_num / len(numeric)
    else:
        out["span_equals_number_in_token"] = None

    out["span_equals_one"] = sum(1 for s in spans if s == 1) / n
    out["span_is_zero"] = sum(1 for s in spans if s == 0) / n
    out["span_exceeds_token_width"] = (
        sum(1 for token, span in records if span > len(token)) / len(records))
    return out


def verdict(rate, tested):
    if tested < MIN_SAMPLES:
        return "UNTESTABLE"
    if rate is None:
        return "UNTESTABLE"
    if rate >= CONFIRM_THRESHOLD:
        return "CONFIRMED"
    if rate <= REFUTE_THRESHOLD:
        return "REFUTED"
    return "WEAK"


# ---------------------------------------------------------------------------
# 阶段 4：外部证据（当数据内无法验证时——去查资料，而不是猜）
# ---------------------------------------------------------------------------
SEARCH_PROMPT = """一个序列数据里出现了这样的记号：{pattern}
（说明：{description}；样例：{example}）

请给出 2-4 条用于检索的查询串，目标是查明：
  1) 什么工具/流程会产生这种记号；
  2) 该工具或相关文献对这个记号的约定是什么。

优先检索"产生它的工具"的官方文档，其次是论文、分析流程文档，再次是论坛/issue。
只输出 JSON：{{"queries": ["...", "..."]}}
"""

CLAIM_PROMPT = """正在判定的记号是：{pattern}（说明：{description}）
样例：{example}

下面是检索返回的结果（每条含 title / url / snippet）。
请判断：这些结果里有没有**明确记载这个记号**的约定？

检索结果：
{hits}

要求：
1. 只依据上面的检索结果，不得引用其中没有的链接；
2. 只接受**与上述记号相关**的内容；若结果讲的是别的记号或别的主题，视为没有记载；
3. source_url 必须是上面出现过的 url，原样照抄；
4. supporting_snippet 必须是上面某条结果里真实出现的文字片段；
5. source_type 从这些值里选一个：official_docs / peer_reviewed / protocol / forum / unknown；
6. 如果没有明确记载**这个记号**，把 claim 留空字符串。

只输出 JSON：
{{"claim": "<一句话结论，或空串>", "source_url": "<上面的某个 url，或空串>",
  "supporting_snippet": "<原文片段，或空串>", "source_type": "<类型>"}}
"""


def _mock_search(_query: str) -> str:
    """离线自检用：返回**明显是占位符**的合成检索结果。

    刻意使用 example.invalid（RFC 2606 保留域名），避免被误当成真实文献。
    """
    return json.dumps([
        {"title": "[合成占位] 某序列工具文档：括号数字记号",
         "url": "https://example.invalid/docs/seq-tool#bracket-markers",
         "snippet": "标记 *[k]* 表示该位置原有 k 个残基未能确定，已在输出中折叠为记号。"},
        {"title": "[合成占位] 论坛讨论：折叠记号的约定",
         "url": "https://example.invalid/forum/thread/12345",
         "snippet": "不同流程对同一记号的处理并不一致，使用时请以产生该文件的工具文档为准。"},
    ], ensure_ascii=False)


def _mock_claim(prompt: str) -> str:
    """离线自检用：只对"检索结果确实记载了该记号"的情况给结论。

    mock 检索结果只记载了 `*[k]*` 记号，因此对其它记号必须返回空
    —— 模拟"查了但没查到"，而不是硬套一条不相关的结论。
    """
    if "star_bracket_number" not in prompt:
        return json.dumps({"claim": "", "source_url": "",
                           "supporting_snippet": "", "source_type": "unknown"},
                          ensure_ascii=False)
    return json.dumps({
        "claim": "该记号表示该处原有 k 个残基未能确定，被折叠为记号",
        "source_url": "https://example.invalid/docs/seq-tool#bracket-markers",
        "supporting_snippet": "标记 *[k]* 表示该位置原有 k 个残基未能确定，已在输出中折叠为记号。",
        "source_type": "official_docs",
    }, ensure_ascii=False)


def call_search(mode, cmd, query):
    if mode == "mock":
        return _mock_search(query)
    if not cmd:
        return None
    try:
        result = subprocess.run(cmd, shell=True, input=query, capture_output=True,
                                text=True, encoding="utf-8", errors="replace")
    except OSError:
        return None
    if result.returncode != 0:
        return None
    return result.stdout


def parse_hits(raw):
    if not raw:
        return []
    text = raw.strip()
    if text.startswith("```"):
        text = re.sub(r"^```[a-zA-Z]*\n?|```$", "", text).strip()
    try:
        data = json.loads(text)
    except json.JSONDecodeError:
        return []
    if not isinstance(data, list):
        return []
    hits = []
    for h in data:
        if isinstance(h, dict) and h.get("url"):
            hits.append({"title": str(h.get("title", "")), "url": str(h["url"]),
                         "snippet": str(h.get("snippet", ""))})
    return hits


def parse_claim(raw, hits):
    """解析模型提取的结论，并**校验它引用的 URL 确实来自检索结果**。

    这是防幻觉的关键一步：模型不能引用一个检索结果里根本不存在的链接。
    """
    text = (raw or "").strip()
    if text.startswith("```"):
        text = re.sub(r"^```[a-zA-Z]*\n?|```$", "", text).strip()
    try:
        claim = json.loads(text)
    except json.JSONDecodeError:
        return None, "模型输出不是合法 JSON"
    if not isinstance(claim, dict):
        return None, "模型输出不是对象"
    # 字段名兼容：机器生成的证据多用 claim.text，模型直出时多用 claim.claim
    text = (claim.get("claim") or claim.get("text") or "").strip()
    if not text:
        return None, "未从检索结果中得到明确结论"
    url = (claim.get("source_url") or "").strip()
    known = {h["url"] for h in hits}
    if url not in known:
        return None, f"引用了检索结果中不存在的来源：{url or '(空)'}"
    snippet = (claim.get("supporting_snippet") or "").strip()
    if not snippet:
        return None, "缺少支撑原文片段"
    if not any(snippet[:20] in h["snippet"] or h["snippet"][:20] in snippet
               for h in hits if h["url"] == url):
        return None, "支撑片段与来源内容不符"
    claim["source_type"] = claim.get("source_type") if claim.get("source_type") in SOURCE_TIERS \
        else "unknown"
    # 归一化字段名，下游统一用 claim
    claim["claim"] = text
    return claim, ""


def external_verdict(claim):
    """外部证据能支撑到什么程度。"""
    if not claim:
        return "UNTESTABLE"
    if claim["source_type"] in DOCUMENTED_SOURCES:
        return "DOCUMENTED"
    return "UNVERIFIED_HYPOTHESIS"


# ---------------------------------------------------------------------------
# 阶段 4 模式 A：待查清单（脚本产出，agent 执行检索）
#
# 为什么不让脚本自己检索：
#   实测把搜索结果页抓下来再解析 HTML 非常脆弱——查询词稍偏，返回的可能是
#   完全不相关的页面（例如搜序列记号却返回视频网站），而且中文搜索引擎对
#   英文技术术语召回很差。真正管用的是：读页面内容、判断相关性、换词重试。
#   这是 agent 的强项，不是脚本的强项。所以脚本只产出"要查什么"，
#   检索由 agent 用它自己的联网工具完成。
# ---------------------------------------------------------------------------
QUERY_MANIFEST_TEMPLATE = {
    "instructions": [
        "对 unresolved 里的每一项，使用你的联网检索能力（web_search / web_fetch 等）查明该记号的约定。",
        "优先找**产生该记号的工具**的官方文档；其次论文；再次分析流程文档；最后论坛/issue。",
        "每条结论必须附带：能打开的 url、该页面上**原文出现**的片段、来源类型。",
        "查不到就如实留空——不要用相似记号的结论套用，不要凭记忆编造来源。",
        "把结果写成 JSON 到 --evidence文件，格式见 response_schema，然后重新运行本脚本并加 --evidence。",
    ],
    "evidence_quality_rules": {
        "DOCUMENTED": "需要 official_docs / peer_reviewed / protocol 级别的来源，"
                      "且 url 与原文片段都要能对应上",
        "UNVERIFIED_HYPOTHESIS": "只有 forum 级来源，或来源无法对应具体结论",
        "UNTESTABLE": "没有找到任何记载，或找到的都是其它记号的记载",
    },
}


def _emit_query_manifest(args, findings, need_external, by_pattern):
    """产出待查清单：告诉 agent 要查什么、什么算合格证据、结果怎么写回来。"""
    unresolved = []
    for f in need_external:
        info = by_pattern[f["pattern"]]
        unresolved.append({
            "pattern": f["pattern"],
            "hypothesis": f.get("hypothesis", ""),
            "description": info["description"],
            "example": info.get("example", ""),
            "occurrences": info.get("count", 0),
            "distinct_tokens": info.get("distinct", []),
            "why_unresolved": f.get("reason") or f["verdict"],
            "question": f"什么工具或流程会产生这种记号？它对“{info['description']}”的约定是什么？",
            "suggested_queries": _suggest_queries(f, info),
        })
    manifest = dict(QUERY_MANIFEST_TEMPLATE)
    manifest.update({
        "sequence_column": args.seq_col,
        "reference_column": args.ref_col,
        "unresolved": unresolved,
        "response_schema": {
            "<pattern>": {
                "search_queries": ["<你实际用过的检索式>"],
                "hits": [{"title": "<页面标题>", "url": "<可打开的链接>",
                          "snippet": "<页面上原文出现的片段>"}],
                "claim": {
                    "text": "<一句话结论；查不到就留空串>",
                    "source_url": "<必须是 hits 里出现过的 url，或空串>",
                    "supporting_snippet": "<必须与 hits 中该 url 的 snippet 相符，或空串>",
                    "source_type": "official_docs|peer_reviewed|protocol|forum|unknown",
                },
            }
        },
    })
    with open(args.emit_queries, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, ensure_ascii=False, indent=2)


def _suggest_queries(finding, info):
    """给出检索式建议。agent 可以改，这只是起点。"""
    desc = info["description"]
    tokens = " ".join(str(t) for t in info.get("distinct", [])[:3])
    return [
        f"{tokens} sequence notation meaning".strip(),
        f"{desc} 序列记号 含义 工具",
        f"which tool produces {tokens} marker in sequence files".strip(),
    ]


# ---------------------------------------------------------------------------
# 阶段 4 模式 A（续）：读入 agent 查回的证据并裁定
# ---------------------------------------------------------------------------
def _ingest_evidence(args, findings, need_external, by_pattern):
    """读入 agent 提供的证据文件，做与 --search-cmd 完全相同的校验与裁定。

    校验逻辑不因来源不同而放松：url 必须在该记号的 hits 里真实出现，
    支撑片段必须与来源相符，否则该结论被拒、降级 UNTESTABLE。
    """
    try:
        # utf-8-sig：容忍 BOM。agent 或编辑器写出的 JSON 常带 BOM，
        # 用纯 utf-8 读会在第一个字符处直接失败。
        with open(args.evidence, encoding="utf-8-sig") as fh:
            evidence = json.load(fh)
    except (OSError, json.JSONDecodeError) as exc:
        print(f"  [错误] 无法读取证据文件：{exc}")
        for f in need_external:
            f["external"] = {"verdict": "UNTESTABLE", "reason": f"证据文件不可读：{exc}"}
        return
    if not isinstance(evidence, dict):
        print("  [错误] 证据文件应为对象：{pattern: {...}}")
        for f in need_external:
            f["external"] = {"verdict": "UNTESTABLE", "reason": "证据文件格式错误"}
        return

    print(f"\n阶段 4  读入 agent 提供的证据：{args.evidence}")
    for f in need_external:
        pattern = f["pattern"]
        entry = evidence.get(pattern)
        if not isinstance(entry, dict):
            f["external"] = {"verdict": "UNTESTABLE", "reason": "证据文件中无该记号的条目"}
            print(f"  · {pattern}: 证据文件中没有对应条目 → UNTESTABLE")
            continue

        hits = [h for h in (entry.get("hits") or []) if isinstance(h, dict) and h.get("url")]
        f["search_queries"] = entry.get("search_queries") or []
        f["search_hits"] = [{"url": h["url"], "title": h.get("title", "")} for h in hits]
        print(f"  · {pattern}: agent 提供 {len(hits)} 条来源")

        claim_raw = entry.get("claim")
        if not isinstance(claim_raw, dict):
            f["external"] = {"verdict": "UNTESTABLE", "reason": "未提供结论对象"}
            print("      未提供结论 → UNTESTABLE")
            continue
        claim, why = parse_claim(json.dumps(claim_raw, ensure_ascii=False), hits)
        ev = external_verdict(claim)
        f["external"] = {
            "verdict": ev, "reason": why,
            "claim": (claim or {}).get("claim", ""),
            "source_url": (claim or {}).get("source_url", ""),
            "source_type": (claim or {}).get("source_type", ""),
            "supporting_snippet": (claim or {}).get("supporting_snippet", ""),
            "provided_by": "agent",
        }
        if claim:
            f["derived_prediction"] = "span_equals_number_in_token"
            print(f"      结论[{ev}] {claim['claim']}")
            print(f"      来源({claim['source_type']}): {claim['source_url']}")
            if ev == "DOCUMENTED":
                print("      注意：文献级证据，未经数据验证；拿到参照列后请转成预测实测。")
        else:
            print(f"      结论被拒（{why}）→ {ev}")


# ---------------------------------------------------------------------------
# 阶段 5：裁定与报告
# ---------------------------------------------------------------------------
def main(argv=None):
    parser = argparse.ArgumentParser(
        description="发现未知数据约定的含义（LLM 提假设，代码做验证）")
    parser.add_argument("csv")
    parser.add_argument("--seq-col", required=True)
    parser.add_argument("--ref-col", default=None,
                        help="未被污染的参照列；没有它则锚定实验无法进行")
    parser.add_argument("--llm", choices=["mock"], default=None,
                        help="用内置 mock 假设离线自检")
    parser.add_argument("--llm-cmd", default=None,
                        help="模型桥命令：从 stdin 读 prompt、回复写 stdout")
    parser.add_argument("--search", choices=["mock"], default=None,
                        help="用内置合成检索结果离线自检（URL 为 example.invalid 占位符）")
    parser.add_argument("--search-cmd", default=None,
                        help="检索桥命令：从 stdin 读查询串、把 JSON 结果列表写 stdout；"
                             "每条为 {title, url, snippet}。适合已有搜索 API 的场景")
    parser.add_argument("--emit-queries", default=None,
                        help="【推荐】产出待查清单 JSON，由具备联网能力的 agent 执行检索，"
                             "结果写成证据文件后用 --evidence 喂回")
    parser.add_argument("--evidence", default=None,
                        help="读入 agent 检索后写好的证据 JSON，做校验与裁定")
    parser.add_argument("--flank", type=int, default=15)
    parser.add_argument("--sample-rows", type=int, default=20000)
    parser.add_argument("--anchor-limit", type=int, default=3000,
                        help="锚定实验最多测多少个记号")
    parser.add_argument("--emit-rules", default=None,
                        help="把确认的结论写成规则建议 JSON")
    args = parser.parse_args(argv)

    if not args.llm and not args.llm_cmd:
        print("[错误] 请指定 --llm mock 或 --llm-cmd '<命令>'")
        return 2
    if args.emit_queries and args.evidence:
        print("[错误] --emit-queries 与 --evidence 不能同时使用："
              "先产出清单、检索后再喂回证据。")
        return 2

    print("阶段 1/5  观察数据（确定性扫描）……")
    profile = observe(args.csv, args.seq_col, args.sample_rows)
    col_profile = profile[args.seq_col]
    print(f"  扫描 {col_profile['rows_scanned']:,} 行")
    print(f"  字符集: {col_profile['alphabet'][:60]}")
    for s in col_profile["suspects"]:
        print(f"  - {s['pattern']}: {s['count']:,} 次  例: {s['example'][:48]}")

    print("\n阶段 2/5  请模型提出可证伪的假设……")
    prompt = PROMPT_TEMPLATE.format(
        predictions="\n".join(f"  - {k}：{v}" for k, v in SUPPORTED_PREDICTIONS.items()),
        profile=json.dumps(col_profile, ensure_ascii=False, indent=1),
    )
    raw = call_llm(args.llm, args.llm_cmd, prompt)
    hypotheses, err = parse_hypotheses(raw)
    if err:
        print(f"  [警告] {err} → 无法提出假设，结束（不猜）")
        return 2
    print(f"  收到 {len(hypotheses)} 条假设")

    print("\n阶段 3/5  跑验证实验（确定性，模型不参与判定）……")
    by_pattern = {s["pattern"]: s for s in col_profile["suspects"]}
    findings = []
    for h in hypotheses:
        pattern = h.get("pattern", "")
        info = by_pattern.get(pattern)
        prediction = h.get("prediction", "")
        entry = {"pattern": pattern, "hypothesis": h.get("hypothesis", ""),
                 "rationale": h.get("rationale", ""), "prediction": prediction}
        if info is None:
            entry.update(verdict="UNTESTABLE", reason="该记号在数据中不存在", tested=0)
            findings.append(entry); continue
        if not prediction:
            entry.update(verdict="UNTESTABLE", reason="模型未给出可验证的预测", tested=0)
            findings.append(entry); continue
        if not args.ref_col:
            entry.update(verdict="UNTESTABLE", reason="未提供参照列，转外部证据阶段", tested=0)
            # 不 continue：先做数据内尝试，失败后由阶段 4 去查资料
            findings.append(entry)
            continue

        probe = probe_anchor_span(args.csv, args.seq_col, args.ref_col,
                                  info["regex"], args.flank, args.anchor_limit)
        rates = evaluate_predictions(probe)
        rate = rates.get(prediction)
        v = verdict(rate, probe["tested"])
        entry.update(
            verdict=v, tested=probe["tested"],
            match_rate=(round(rate, 4) if rate is not None else None),
            span_examples=probe["records"][:3],
            all_rates={k: (round(x, 4) if x is not None else None)
                       for k, x in rates.items()},
        )
        findings.append(entry)

    # ---------------- 阶段 4：外部证据（数据内判不了的，去查资料，而不是猜） -------------
    #
    # 三种模式，按推荐程度排序：
    #   A. --emit-queries FILE  → 脚本只产出"待查清单"，由 agent（具备 web_search/web_fetch
    #      的宿主，如 DSH）去检索并写回证据，再 --evidence FILE 回来裁定。
    #      这是推荐用法：检索需要"读页面、判断相关性、换词重试"，正是 agent 的强项，
    #      而脚本拼 URL 抓搜索结果页再解析 HTML 非常脆弱（实测会返回完全不相关的结果）。
    #   B. --search mock       → 离线合成检索，仅用于跑通与单测。
    #   C. --search-cmd "..."  → 自带检索桥，全自动；适合已有搜索引擎 API 的场景。
    need_external = [f for f in findings
                     if f["verdict"] in ("UNTESTABLE", "WEAK")
                     and f.get("pattern") in by_pattern]

    if need_external and args.emit_queries:
        _emit_query_manifest(args, findings, need_external, by_pattern)
        print(f"\n阶段 4  已产出待查清单：{args.emit_queries}")
        print("  下一步：由具备联网能力的 agent 按清单检索，把结果写成证据文件，")
        print("  再执行：python discover_conventions.py <同一 csv> --seq-col ... "
              f"--evidence <证据文件> --emit-rules rules.json")
        print("  说明：检索需要读页面、判断相关性、换词重试，这是 agent 的职责；")
        print("        脚本只负责产出清单与最终裁定，不做网络请求。")
        return 0

    if need_external and args.evidence:
        _ingest_evidence(args, findings, need_external, by_pattern)
    elif need_external and not (args.search or args.search_cmd):
        print(f"\n阶段 4  数据内无法判定 {len(need_external)} 条")
        print("  未提供外部证据来源，这些保持 UNTESTABLE —— 不猜。")
        print("  建议：加 --emit-queries q.json 产出待查清单，让 agent 去查。")
    elif need_external:
        print(f"\n阶段 4  数据内无法判定 {len(need_external)} 条，转入外部证据检索……")
        for f in need_external:
            info = by_pattern[f["pattern"]]
            q_prompt = SEARCH_PROMPT.format(
                pattern=f["pattern"], description=info["description"],
                example=info.get("example", "")[:80])
            raw_q = call_llm(args.llm, args.llm_cmd, q_prompt)
            try:
                queries = json.loads(re.sub(r"^```[a-zA-Z]*\n?|```$", "", raw_q.strip())).get("queries", [])
            except (json.JSONDecodeError, AttributeError):
                queries = []
            queries = [q for q in queries if isinstance(q, str) and q.strip()][:4]
            f["search_queries"] = queries
            if not queries:
                f["external"] = {"verdict": "UNTESTABLE", "reason": "模型未给出检索式"}
                print(f"  · {f['pattern']}: 模型未给出检索式 → UNTESTABLE")
                continue

            hits = []
            for q in queries:
                hits.extend(parse_hits(call_search(args.search, args.search_cmd, q)))
            # 去重
            seen, uniq = set(), []
            for h in hits:
                if h["url"] not in seen:
                    seen.add(h["url"]); uniq.append(h)
            hits = uniq
            f["search_hits"] = [{"url": h["url"], "title": h["title"]} for h in hits]
            print(f"  · {f['pattern']}: 检索式 {len(queries)} 条 → 命中 {len(hits)} 条")
            if not hits:
                f["external"] = {"verdict": "UNTESTABLE", "reason": "检索无命中"}
                print("      无命中 → UNTESTABLE（不猜）")
                continue

            claim_raw = call_llm(args.llm, args.llm_cmd, CLAIM_PROMPT.format(
                pattern=f["pattern"], description=info["description"],
                example=info.get("example", "")[:80],
                hits=json.dumps(hits, ensure_ascii=False, indent=1)))
            claim, why = parse_claim(claim_raw, hits)
            ev = external_verdict(claim)
            f["external"] = {"verdict": ev, "reason": why,
                             "claim": (claim or {}).get("claim", ""),
                             "source_url": (claim or {}).get("source_url", ""),
                             "source_type": (claim or {}).get("source_type", ""),
                             "supporting_snippet": (claim or {}).get("supporting_snippet", "")}
            if claim:
                # 外部证据若指向"计数"，就把它转成一个**可检验的预测**，等日后有参照列时验证
                f["derived_prediction"] = "span_equals_number_in_token"
                print(f"      结论[{ev}] {claim['claim']}")
                print(f"      来源({claim['source_type']}): {claim['source_url']}")
                if ev == "DOCUMENTED":
                    print("      注意：这是文献级证据，未经数据验证；"
                          "一旦拿到参照列请转成预测实测。")
            else:
                print(f"      未能提取可引用结论（{why}）→ {ev}")

    print("\n阶段 5  裁定")
    print("=" * 104)
    print(f"{'记号':<22}{'假设':<12}{'数据内判定':<13}{'样本':>7}{'命中率':>8}   "
          f"{'外部证据':<22}{'最终':>12}")
    print("-" * 104)
    for f in findings:
        rate = f.get("match_rate")
        ext = f.get("external", {})
        ext_s = f"{ext.get('verdict','')}"
        if ext.get("source_type"):
            ext_s += f"({ext['source_type']})"
        final = _final_verdict(f)
        print(f"{f['pattern']:<22}{str(f.get('hypothesis',''))[:10]:<12}"
              f"{f['verdict']:<13}{f.get('tested',0):>7}"
              f"{(f'{100*rate:.1f}%' if rate is not None else '—'):>8}   "
              f"{ext_s[:20]:<22}{final:>12}")
    print("-" * 104)
    for key, meaning in VERDICT_MEANING.items():
        print(f"  {key:<22}{meaning}")
    print()
    for f in findings:
        ext = f.get("external", {})
        if ext.get("source_url"):
            print(f"  · {f['pattern']} 外部来源：{ext['source_url']}")
            print(f"      「{ext.get('supporting_snippet','')[:90]}」")
        if f.get("reason"):
            print(f"  · {f['pattern']}/{f.get('hypothesis','')}: {f['verdict']} —— {f['reason']}")
    for f in findings:
        if f.get("search_queries"):
            print(f"  · {f['pattern']} 检索式：{f['search_queries']}")

    if args.emit_rules:
        rules = _suggest_rules(findings)
        with open(args.emit_rules, "w", encoding="utf-8") as fh:
            json.dump({"source": args.csv, "sequence_column": args.seq_col,
                       "reference_column": args.ref_col,
                       "findings": findings, "suggested_rules": rules,
                       "note": "建议需人工确认后再交给 curate.py。"
                               "CONFIRMED 可较放心采用；DOCUMENTED 为文献级证据，"
                               "建议先在小样本上验证；UNVERIFIED_HYPOTHESIS 与 UNTESTABLE 不得自动生效。"},
                      fh, ensure_ascii=False, indent=2)
        print(f"\n规则建议已写入: {args.emit_rules}（需人工确认）")
    return 0


def _final_verdict(f):
    """数据内证据优先于外部证据。"""
    if f["verdict"] == "CONFIRMED":
        return "CONFIRMED"
    if f["verdict"] == "REFUTED":
        return "REFUTED"
    ext = f.get("external", {}).get("verdict")
    if ext in ("DOCUMENTED", "UNVERIFIED_HYPOTHESIS"):
        return ext
    return f["verdict"]


def _suggest_rules(findings):
    """把已确认的结论翻译成 curate.py 可选用的规则建议。

    只有 CONFIRMED（数据内验证）与 DOCUMENTED（文献级）会进入建议；
    后者标注为需先在小样本上复核。
    """
    rules = []
    for f in findings:
        final = _final_verdict(f)
        if final not in ("CONFIRMED", "DOCUMENTED"):
            continue
        pred = f.get("prediction") or f.get("derived_prediction")
        evidence = (f"数据内锚定 {f['tested']} 例，命中率 {f['match_rate']}"
                    if final == "CONFIRMED" else
                    f"文献级证据：{f.get('external', {}).get('source_url', '')}")
        needs_check = final == "DOCUMENTED"
        if f["pattern"] == "star_bracket_number" and pred == "span_equals_number_in_token":
            rules.append({"regex": r"\*\[ *(\d+)\]\*", "action": "expand_to_n_by_number",
                          "evidence_tier": final, "evidence": evidence,
                          "requires_small_sample_check": needs_check,
                          "maps_to": "占位符 -> k 个 N"})
        elif pred == "span_equals_token_width" and f["pattern"] == "lowercase_letters":
            rules.append({"regex": "[a-z]", "action": "uppercase_keep",
                          "evidence_tier": final, "evidence": evidence,
                          "requires_small_sample_check": needs_check,
                          "maps_to": "小写 -> 同一碱基的大写（非插入记号）"})
    return rules


if __name__ == "__main__":
    sys.exit(main())
