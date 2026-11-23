"""发现循环回归测试：观察 → LLM 提假设 → 代码验证 → 裁定。

重点验证这个模块的**纪律**，而不只是它能跑：
  - 模型只能从受支持的预测清单里选；
  - 模型给出不受支持的预测 → 置空 → 该假设判 UNTESTABLE，绝不当成通过；
  - 样本不足 → UNTESTABLE（不允许小样本下结论）；
  - 参照列缺失 → UNTESTABLE（不给"看起来像"留后门）；
  - 侧翼只收字母：混入噪声时不得把样本数静默变成 0（这是本项目踩过的真实缺陷）。

全部离线运行，不连任何模型。
"""
import contextlib
import io
import json
import os
import sys

from discover_conventions import (
    DOCUMENTED_SOURCES,
    MIN_SAMPLES,
    SUPPORTED_PREDICTIONS,
    VERDICT_MEANING,
    evaluate_predictions,
    external_verdict,
    main,
    observe,
    parse_claim,
    parse_hypotheses,
    probe_anchor_span,
    verdict,
)

HERE = os.path.dirname(os.path.abspath(__file__))
DEMO = os.path.join(HERE, "examples", "convention_discovery_demo.csv")
ANCHOR = os.path.join(HERE, "examples", "placeholder_anchor_demo.csv")

FAILURES = []


def check(name, condition):
    print(("PASS  " if condition else "FAIL  ") + name)
    if not condition:
        FAILURES.append(name)


def run_cli(argv):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        code = main(argv)
    return code, buf.getvalue()


def table_part(out):
    """只取裁定表格的表格体（两条长分隔线之间）。

    不能用"图例文字"来切分：等级说明里本来就写着 CONFIRMED 等字样，
    切不干净会把图例当成表格行，导致断言失真。
    """
    lines = out.splitlines()
    dashes = [i for i, line in enumerate(lines)
              if line.strip() and set(line.strip()) == {"-"} and len(line.strip()) > 20]
    if len(dashes) >= 2:
        return "\n".join(lines[dashes[0] + 1:dashes[1]])
    return out


# ---------------------------------------------------------------- 受支持的预测是封闭集合
check("受支持的预测清单非空且命名稳定",
      "span_equals_number_in_token" in SUPPORTED_PREDICTIONS
      and "span_equals_token_width" in SUPPORTED_PREDICTIONS)

# ---------------------------------------------------------------- 输出解析的纪律
ok, err = parse_hypotheses(json.dumps([
    {"pattern": "star_bracket_number", "hypothesis": "计数",
     "prediction": "span_equals_number_in_token"}]))
check("合法 JSON 可解析", len(ok) == 1 and not err
      and ok[0]["prediction"] == "span_equals_number_in_token")

bad, err2 = parse_hypotheses(json.dumps([
    {"pattern": "star_bracket_number", "hypothesis": "计数",
     "prediction": "我用眼睛看出来的"}]))
check("不受支持的预测名被置空并留痕",
      len(bad) == 1 and bad[0]["prediction"] == ""
      and "不受支持" in bad[0]["rationale"])

not_json, err3 = parse_hypotheses("我觉得是计数")
check("非 JSON 输出被拒绝且给出原因",
      not_json == [] and "不是合法 JSON" in err3)

# ---------------------------------------------------------------- 裁定阈值
check("样本不足不得下结论", verdict(1.0, MIN_SAMPLES - 1) == "UNTESTABLE")
check("高命中率判 CONFIRMED", verdict(1.0, 1000) == "CONFIRMED")
check("低命中率判 REFUTED", verdict(0.0, 1000) == "REFUTED")
check("命中率缺失判 UNTESTABLE", verdict(None, 1000) == "UNTESTABLE")

# ---------------------------------------------------------------- 观察阶段
profile = observe(DEMO, "seq_with_markers", 1000)
suspects = {s["pattern"]: s for s in profile["seq_with_markers"]["suspects"]}
check("观察到 *[k]* 记号", "star_bracket_number" in suspects
      and suspects["star_bracket_number"]["count"] > 0)
check("观察到小写字母", "lowercase_letters" in suspects
      and suspects["lowercase_letters"]["count"] > 0)

# ---------------------------------------------------------------- 侧翼只收字母（真实缺陷回归）
seq = "AAAA*[3]*GGGG"
probe = probe_anchor_span(DEMO, "seq_with_markers", "reference_sequence",
                          suspects["lowercase_letters"]["regex"], 15, 500)
check("侧翼不受记号污染：小写样本数 > 0（曾静默为 0）", probe["tested"] > 0)

# ---------------------------------------------------------------- 端到端（合成演示数据）
code, out = run_cli([DEMO, "--seq-col", "seq_with_markers",
                     "--ref-col", "reference_sequence", "--llm", "mock"])
check("发现循环端到端返回 0", code == 0)
check("计数假设被 CONFIRMED", "CONFIRMED" in table_part(out))
check("宽度假设未被误判为成立（应 WEAK 或 REFUTED）",
      "WEAK" in table_part(out) or "REFUTED" in table_part(out))
check("小写假设被 CONFIRMED 且样本非零",
      "lowercase_letters" in out and "0        —" not in out)

# ------------------------------------------- 缺参照列且未接检索：一律 UNTESTABLE
code2, out2 = run_cli([DEMO, "--seq-col", "seq_with_markers", "--llm", "mock"])
check("无参照列且未接检索接口时全部 UNTESTABLE（不给猜测留后门）",
      code2 == 0 and "UNTESTABLE" in table_part(out2)
      and "CONFIRMED" not in table_part(out2)
      and "REFUTED" not in table_part(out2)
      and "DOCUMENTED" not in table_part(out2))

# ---------------------------------------------------------------- 样本太少 → 不裁定
code3, out3 = run_cli([ANCHOR, "--seq-col", "seq_with_placeholder",
                       "--ref-col", "reference_sequence", "--llm", "mock"])
check("仅 3 行样本时不下结论",
      code3 == 0 and "UNTESTABLE" in table_part(out3)
      and "CONFIRMED" not in table_part(out3))

# ---------------------------------------------------------------- 外部证据层
HITS = [
    {"title": "工具文档", "url": "https://real.example/docs#marker",
     "snippet": "标记 *[k]* 表示该位置原有 k 个残基未能确定，已折叠为记号。"},
    {"title": "论坛贴", "url": "https://real.example/forum/1",
     "snippet": "这个记号在不同流程里含义不一致，请以工具文档为准。"},
]

# 防幻觉核心：模型不得引用检索结果里不存在的链接
ghost = json.dumps({"claim": "它是计数", "source_url": "https://made-up.example/x",
                    "supporting_snippet": "标记 *[k]* 表示该位置原有 k 个残基未能确定，已折叠为记号。",
                    "source_type": "official_docs"}, ensure_ascii=False)
claim, why = parse_claim(ghost, HITS)
check("引用不存在的来源被拒绝（防编造引用）", claim is None and "不存在" in why)

# 支撑片段必须真的来自来源
mismatch = json.dumps({"claim": "它是计数", "source_url": "https://real.example/docs#marker",
                       "supporting_snippet": "文档里根本没写这句话。",
                       "source_type": "official_docs"}, ensure_ascii=False)
claim2, why2 = parse_claim(mismatch, HITS)
check("支撑片段与来源不符被拒绝", claim2 is None and "不符" in why2)

# 合法引用被接受，并按来源类型定级
good = json.dumps({"claim": "k 是被折叠的残基数",
                   "source_url": "https://real.example/docs#marker",
                   "supporting_snippet": "标记 *[k]* 表示该位置原有 k 个残基未能确定，已折叠为记号。",
                   "source_type": "official_docs"}, ensure_ascii=False)
claim3, why3 = parse_claim(good, HITS)
check("合法引用被接受", claim3 is not None and not why3)
check("官方文档 → DOCUMENTED", external_verdict(claim3) == "DOCUMENTED")

forum = dict(claim3, source_type="forum")
check("论坛来源降级为 UNVERIFIED_HYPOTHESIS（可靠性不足）",
      external_verdict(forum) == "UNVERIFIED_HYPOTHESIS")
check("未提取到结论 → UNTESTABLE", external_verdict(None) == "UNTESTABLE")
check("DOCUMENTED 与 CONFIRMED 在清单中是不同等级",
      "DOCUMENTED" in VERDICT_MEANING and "CONFIRMED" in VERDICT_MEANING
      and "official_docs" in DOCUMENTED_SOURCES)

# 端到端：无参照列时走外部证据，且不同记号得到不同结论
code4, out4 = run_cli([DEMO, "--seq-col", "seq_with_markers",
                       "--llm", "mock", "--search", "mock"])
check("无参照列时仍能给出结论（不再一律 UNTESTABLE）",
      code4 == 0 and "DOCUMENTED" in table_part(out4))
check("未查到记载的记号仍为 UNTESTABLE（不硬套）",
      "UNTESTABLE" in table_part(out4))
check("外部来源与其原文片段被打印出来",
      "real" not in out4 and "example.invalid" in out4 and "残基未能确定" in out4)

# ------------------------------------------- 两阶段 agent 模式：清单 → 证据 → 裁定
import tempfile

TMP = tempfile.mkdtemp(prefix="discover_test_")
QFILE = os.path.join(TMP, "queries.json")
EFILE = os.path.join(TMP, "evidence.json")

# 第一趟：产出待查清单
code5, out5 = run_cli([DEMO, "--seq-col", "seq_with_markers", "--llm", "mock",
                       "--emit-queries", QFILE])
check("待查清单模式返回 0", code5 == 0)
check("清单文件已生成", os.path.isfile(QFILE))
with open(QFILE, encoding="utf-8-sig") as f:
    manifest = json.load(f)
check("清单含未解决记号与检索式建议",
      manifest["unresolved"] and manifest["unresolved"][0]["suggested_queries"])
check("检索式建议含代码搜索优先顺序（实测网页搜索对记号无效）",
      "code_search" in manifest["unresolved"][0]["suggested_queries"]
      and "search_strategy" in manifest
      and manifest["search_strategy"]["first_choice"].startswith("代码搜索"))
check("清单含证据质量规则与回填 schema",
      "DOCUMENTED" in manifest["evidence_quality_rules"]
      and "response_schema" in manifest)
check("清单说明由 agent 执行检索",
      any("联网" in s for s in manifest["instructions"]))

# 第二趟：agent 写回证据（含 BOM，模拟真实写文件）
evidence = {
    "star_bracket_number": {
        "search_queries": ["marker notation"],
        "hits": [{"title": "工具文档", "url": "https://real.example/docs#m",
                  "snippet": "标记 *[k]* 表示该位置原有 k 个残基未能确定，已折叠。"}],
        # 故意用 text 字段（文档 schema 的写法），验证字段名兼容
        "claim": {"text": "k 是被折叠的残基数", "source_url": "https://real.example/docs#m",
                  "supporting_snippet": "标记 *[k]* 表示该位置原有 k 个残基未能确定，已折叠。",
                  "source_type": "official_docs"},
    },
    "lowercase_letters": {
        "search_queries": ["lowercase"],
        "hits": [{"title": "无关页", "url": "https://real.example/other",
                  "snippet": "与本记号无关的内容。"}],
        "claim": {"text": "小写代表低置信度", "source_url": "https://real.example/other",
                  "supporting_snippet": "小写代表低置信度", "source_type": "official_docs"},
    },
}
with open(EFILE, "w", encoding="utf-8-sig") as f:
    json.dump(evidence, f, ensure_ascii=False)

code6, out6 = run_cli([DEMO, "--seq-col", "seq_with_markers", "--llm", "mock",
                       "--evidence", EFILE])
body6 = table_part(out6)
check("证据喂回模式返回 0", code6 == 0)
check("带 BOM 的证据文件可读（真实写文件常见）", "无法读取" not in out6)
check("claim.text 字段名被兼容", "DOCUMENTED" in body6)
check("片段与来源不符的结论被拒", body6.count("UNTESTABLE") >= 1)
check("引用的真实来源被打印", "real.example/docs#m" in out6)

# 互斥参数
code7, out7 = run_cli([DEMO, "--seq-col", "seq_with_markers", "--llm", "mock",
                       "--emit-queries", QFILE, "--evidence", EFILE])
check("--emit-queries 与 --evidence 互斥报错", code7 == 2)

print()
if FAILURES:
    raise SystemExit(f"失败 {len(FAILURES)} 项: {FAILURES}")
print("所有断言通过")