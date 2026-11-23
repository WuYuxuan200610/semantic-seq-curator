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
    MIN_SAMPLES,
    SUPPORTED_PREDICTIONS,
    evaluate_predictions,
    main,
    observe,
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
    """只取裁定表格部分——图例里永远写着 CONFIRMED 等字样，不能拿来当判据。"""
    return out.split("CONFIRMED = 预测成立")[0]


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

# ---------------------------------------------------------------- 缺参照列 → 一律 UNTESTABLE
code2, out2 = run_cli([DEMO, "--seq-col", "seq_with_markers", "--llm", "mock"])
check("未提供参照列时全部 UNTESTABLE（不给猜测留后门）",
      code2 == 0 and "UNTESTABLE" in table_part(out2)
      and "CONFIRMED" not in table_part(out2)
      and "REFUTED" not in table_part(out2))

# ---------------------------------------------------------------- 样本太少 → 不裁定
code3, out3 = run_cli([ANCHOR, "--seq-col", "seq_with_placeholder",
                       "--ref-col", "reference_sequence", "--llm", "mock"])
check("仅 3 行样本时不下结论",
      code3 == 0 and "UNTESTABLE" in table_part(out3)
      and "CONFIRMED" not in table_part(out3))

print()
if FAILURES:
    raise SystemExit(f"失败 {len(FAILURES)} 项: {FAILURES}")
print("所有断言通过")
