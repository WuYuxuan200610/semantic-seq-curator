"""目标导向对比：四种做法产出的数据，能不能进模型输入口、还剩多少可用数据。

为什么用这两个指标，而不是"保留了多少信息"：
    训练语言模型唯一的硬要求是——输出里每个字符都在 tokenizer 词表里。
    一个字符喂不进去，这条数据就是坏的。所以第一指标是**合法率**。
    但"全删光"也能拿到 100% 合法率，所以第二指标是**数据产出**（还剩多少字符可用）。
    两者必须同时看。

对比四种真实规则（各自的实现就在本文件里，可逐行核对）：
    ours      本仓库 curate.py 的规则：小写→大写；gap/模糊码→N；占位符 *[k]*→k 个 N
    cleanfa   seismic cleanfa 的规则：模糊码→N；非 IUPAC 字符报错/丢弃（gap、标记字符会丢）
    dumb      最朴素的脚本：只保留 ACGTU，其余全删
    snip      SNIP 的规则：白名单只有 20 种标准氨基酸，含其它字符就整条删除

用法：
    python benchmark_model_ready.py <raw.csv> --seq-col Sequence
    python benchmark_model_ready.py <raw.csv> --seq-col Sequence \
        --out benchmark/model_ready.json --chart benchmark/model_ready.png
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys

PLACEHOLDER_RE = re.compile(r"\*\[ *(\d+)\]\*")
AMBIG = set("NRYWSKMBDHVX")
STD = set("ACGU")
AA20 = set("ARNDCQEGHILKMFPSTWYV")


def rule_ours(seq, alphabet):
    """curate.py 的模型就绪规则。"""
    out, i = [], 0
    while i < len(seq):
        m = PLACEHOLDER_RE.match(seq, i)
        if m:
            out.append("N" * int(m.group(1)))
            i = m.end()
            continue
        up = seq[i].upper()
        out.append(up if up in STD else "N")
        i += 1
    return "".join(out)


def rule_cleanfa(seq, alphabet):
    """seismic cleanfa 的规则。"""
    out = []
    for ch in seq:
        up = ch.upper()
        if up in STD:
            out.append(up)
        elif up in AMBIG:
            out.append("N")
        # 其它字符（gap、占位符符号、空格）不被接受 → 丢弃
    return "".join(out)


def rule_dumb(seq, alphabet):
    """最朴素的删除式脚本。"""
    return "".join(ch.upper() for ch in seq if ch.upper() in STD)


def main(argv=None):
    parser = argparse.ArgumentParser(description="合法率 / 数据产出 对比")
    parser.add_argument("csv", help="原始 CSV")
    parser.add_argument("--seq-col", default=None, help="序列列（默认自动挑最脏的一列）")
    parser.add_argument("--alphabet", default="ACGUN",
                        help="目标字母表，默认 ACGUN；用于判定合法率")
    parser.add_argument("--out", default=None, help="把结果写成 JSON")
    parser.add_argument("--chart", default=None, help="顺便画一张对比图（需 matplotlib）")
    parser.add_argument("--label", default="本地原始数据",
                        help="图上标注的数据来源说明（不含任何具体数据名）")
    args = parser.parse_args(argv)

    alphabet = set(args.alphabet)
    try:
        handle = open(args.csv, encoding="utf-8-sig", errors="replace", newline="")
    except OSError as exc:
        print(f"[错误] 无法读取：{exc}")
        return 2

    with handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        seq_col = args.seq_col
        if seq_col is None:
            best, best_score = None, -1
            for col in fieldnames:
                joined = ""
                handle.seek(0)
                probe = csv.DictReader(handle)
                for i, row in enumerate(probe):
                    if i >= 300:
                        break
                    joined += row.get(col) or ""
                if len(joined) < 200:
                    continue
                score = sum(1 for ch in joined if ch.upper() not in STD and ch != "-")
                if score > best_score:
                    best, best_score = col, score
            seq_col = best
        if not seq_col or seq_col not in fieldnames:
            print(f"[错误] 未找到序列列，请用 --seq-col 指定。可用：{fieldnames[:15]}……")
            return 2
        handle.seek(0)
        reader = csv.DictReader(handle)
        stats = {
            "本仓库规则（模型就绪）": dict(valid=0, total=0, orig=0, kept=0, rows=0),
            "seismic cleanfa": dict(valid=0, total=0, orig=0, kept=0, rows=0),
            "纯 Python 哑程序": dict(valid=0, total=0, orig=0, kept=0, rows=0),
            "SNIP（蛋白白名单）": dict(valid=0, total=0, orig=0, kept=0, rows=0),
        }
        for row in reader:
            seq = row.get(seq_col) or ""
            if not seq:
                continue
            for name, rule in (("本仓库规则（模型就绪）", rule_ours),
                               ("seismic cleanfa", rule_cleanfa),
                               ("纯 Python 哑程序", rule_dumb)):
                out = rule(seq, alphabet)
                r = stats[name]
                r["rows"] += 1
                r["orig"] += len(seq)
                r["total"] += len(out)
                r["valid"] += sum(1 for c in out if c in alphabet)
                if out:
                    r["kept"] += 1
            r = stats["SNIP（蛋白白名单）"]
            r["rows"] += 1
            r["orig"] += len(seq)
            if all(ch.upper() in AA20 for ch in seq):
                r["kept"] += 1
                r["total"] += len(seq)
                r["valid"] += sum(1 for ch in seq.upper() if ch in alphabet)

    best_yield = max((100 * r["total"] / max(1, r["orig"])) for r in stats.values())
    results = {}
    print("=" * 90)
    print("目标导向对比：输出能不能进模型输入口 / 还剩多少可用数据")
    print("=" * 90)
    print(f"序列列 {seq_col} / 目标字母表 {''.join(sorted(alphabet))}")
    print()
    print(f"{'方法':<24}{'合法率':>10}{'输出字符数':>14}{'净变化':>14}{'序列保留':>14}")
    print("-" * 90)
    for name, r in stats.items():
        if r["total"] == 0:
            valid_s = "—"
        else:
            valid_s = f"{100 * r['valid'] / r['total']:.2f}%"
        fully_valid = r["total"] > 0 and r["valid"] == r["total"]
        top = r["total"] > 0 and abs(100 * r["total"] / r["orig"] - best_yield) < 0.05
        note = ""
        if fully_valid and top:
            note = "  ← 既 100% 合法、又留下最多数据"
        elif r["total"] == 0:
            note = "  ← 整条丢弃，无可用数据"
        net = r["total"] - r["orig"]
        net_s = f"{net:+,}" if r["total"] else "—"
        total_s = f"{r['total']:,}" if r["total"] else "—"
        kept_s = f"{r['kept']:,}/{r['rows']:,}"
        print(f"{name:<24}{valid_s:>10}{total_s:>14}{net_s:>14}{kept_s:>14}{note}")
        results[name] = {
            "valid_rate_pct": round(100 * r["valid"] / r["total"], 2) if r["total"] else None,
            "output_chars": r["total"],
            "input_chars": r["orig"],
            "net_chars": net,
            "sequences_kept": r["kept"], "sequences_total": r["rows"],
            "usable_chars": r["valid"],
        }
    print("-" * 90)
    print("合法率     = 输出字符 ∈ 目标字母表 的比例（<100% 就有字符喂不进模型）")
    print("净变化     = 输出字符数 - 输入字符数")
    print("             >0 不是错误：输入本身已被上游删掉了一部分碱基（占位符就是痕迹），")
    print("             还原后比输入更长；而 cleanfa / 哑程序会把这段长度丢掉，净变化为负。")
    print("序列保留   = 未被整条丢弃的序列数")

    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as dst:
            json.dump({"sequence_column": seq_col,
                       "alphabet": "".join(sorted(alphabet)),
                       "results": results}, dst, ensure_ascii=False, indent=2)
        print(f"\n结果已写入: {args.out}")

    if args.chart:
        _draw(results, args.chart, args.label)
    return 0


def _draw(results, path, label):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False

    names = ["本仓库规则\n(模型就绪)", "seismic\ncleanfa", "纯 Python\n哑程序", "SNIP\n(蛋白白名单)"]
    keys = ["本仓库规则（模型就绪）", "seismic cleanfa", "纯 Python 哑程序", "SNIP（蛋白白名单）"]
    valid = [results[k]["valid_rate_pct"] or 0.0 for k in keys]
    # 图上直接用净变化字符数（比 >100% 的百分比直观，不会被误读成注水）
    net = [results[k]["net_chars"] for k in keys]

    x = np.arange(len(names))
    w = 0.38
    fig, ax = plt.subplots(figsize=(10.5, 5.2), dpi=150)
    ax2 = ax.twinx()
    b1 = ax.bar(x - w / 2, valid, w, label="输出合法率（%，左轴）", color="#2171b5")
    b2 = ax2.bar(x + w / 2, net, w, label="字符净变化（相对输入，右轴）", color="#9ecae1")
    for bar in b1:
        v = bar.get_height()
        ax.annotate(f"{v:.1f}%", (bar.get_x() + bar.get_width() / 2, v),
                    xytext=(0, 3), textcoords="offset points",
                    ha="center", va="bottom", fontsize=9)
    for bar in b2:
        v = bar.get_height()
        ax2.annotate(f"{v:+,}", (bar.get_x() + bar.get_width() / 2, v),
                     xytext=(0, 3 if v >= 0 else -12), textcoords="offset points",
                     ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("输出合法率（%）")
    ax2.set_ylabel("字符净变化")
    ax.set_ylim(0, 118)
    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=9)
    ax.set_title(f"清理输出能否直接训练（基于{label}，未随仓库提供）", fontsize=12)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="lower left", fontsize=9, framealpha=0.95)
    ax.grid(axis="y", linestyle="--", alpha=0.35)
    ax.spines["top"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, bbox_inches="tight")
    print(f"图已写入: {path}")


if __name__ == "__main__":
    sys.exit(main())
