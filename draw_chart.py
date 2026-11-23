"""绘制对比评测柱状图（benchmark/comparison.png）。

读取 benchmark/results.json，生成"方法 × 指标"分组柱状图，供 GitHub README 引用。
"""

import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

HERE = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(HERE, "benchmark", "results.json"), encoding="utf-8") as f:
    results = json.load(f)

# 固定展示顺序（本 Agent 放最后并高亮）
order = ["白板 Agent", "简单 Prompt", "seismic cleanfa", "SNIP", "本 Agent"]
methods = [m for m in order if m in results]
metrics = [
    ("retention", "信息保留率"),
    ("noise_removal", "噪声清除率"),
    ("explainability", "可解释率"),
    ("accuracy", "端到端正确率"),
]

import numpy as np
x = np.arange(len(methods))
n = len(metrics)
width = 0.8 / n
colors = ["#c7c7c7", "#9ecae1", "#fdae6b", "#a1d99b", "#2171b5"]
# 把"本 Agent"单独用深蓝高亮
bar_colors = ["#b0b0b0", "#9ecae1", "#fdae6b", "#a1d99b", "#2171b5"]

fig, ax = plt.subplots(figsize=(11, 5.5), dpi=150)
for i, (key, label) in enumerate(metrics):
    vals = [results[m][key] * 100 for m in methods]
    bars = ax.bar(x + (i - (n - 1) / 2) * width, vals, width, label=label,
                  color=[bar_colors[j] for j in range(len(methods))])
    for bar, v in zip(bars, vals):
        ax.annotate(f"{v:.0f}", (bar.get_x() + bar.get_width() / 2, v),
                    xytext=(0, 3), textcoords="offset points", ha="center", va="bottom",
                    fontsize=8)

ax.set_ylabel("得分（%）")
ax.set_xticks(x)
ax.set_xticklabels(methods, fontsize=10)
ax.set_ylim(0, 115)
ax.set_title("生物序列非标准字符清理：本 Agent vs 基线 / 竞品（金标准 11 字符）", fontsize=12)
ax.legend(loc="upper left", fontsize=9, ncol=2, framealpha=0.9)
ax.grid(axis="y", linestyle="--", alpha=0.35)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

out = os.path.join(HERE, "benchmark", "comparison.png")
plt.tight_layout()
plt.savefig(out, bbox_inches="tight")
print("saved:", out)
