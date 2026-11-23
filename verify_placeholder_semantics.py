"""占位符 *[k]* 语义的实证脚本 —— 本仓库唯一一个"别处可能没有"的方法。

问题：上游工具把碱基换成了 *[k]*。k 到底是"丢失了 k 个碱基"（计数），
      还是"第 k 号标记"（编号）？两者的处理方式完全不同：
        计数 → 还原成 k 个 N，长度回到真实值，可训练；
        编号 → 标记只是一个注解，长度不该变。
      光看数据本身看不出来，因为结构注释列也是从**已损坏**的序列算出来的——
      那是循环证据（本仓库在这个坑上栽过一次）。

方法：找一个**未被污染**的参照列（例如原始 DNA 提取序列），取占位符两侧各 N bp 作锚点，
      在参照列里定位两侧锚点，量出它们之间的真实跨度。这就能直接判定 k 的含义。

用法：
    python verify_placeholder_semantics.py <raw.csv> \
        --seq-col Sequence --ref-col FASTA_sequence [--flank 15]

输出：每个占位符的 (k, 标记宽度, 真实跨度) 抽样，以及
      "真实跨度 == k" 与 "真实跨度 == 标记宽度" 的比例。
"""
from __future__ import annotations

import argparse
import csv
import re
import sys

PLACEHOLDER_RE = re.compile(r"\*\[ *(\d+)\]\*")


def rna_to_dna(text: str) -> str:
    return text.replace("U", "T").replace("u", "t")


def main(argv=None):
    parser = argparse.ArgumentParser(description="实证占位符 *[k]* 中 k 的含义")
    parser.add_argument("csv", help="原始 CSV 路径")
    parser.add_argument("--seq-col", default="Sequence", help="含占位符的序列列")
    parser.add_argument("--ref-col", default="FASTA_sequence",
                        help="未被污染的参照序列列（用于侧翼锚定）")
    parser.add_argument("--flank", type=int, default=15,
                        help="侧翼锚点长度，默认 15；太短会在参照里误匹配")
    parser.add_argument("--show", type=int, default=10, help="抽样打印条数")
    args = parser.parse_args(argv)

    if args.flank < 8:
        print("[错误] --flank 至少 8；锚点太短会误匹配。")
        return 2

    match_k = match_width = tested = sampled = 0
    samples = []

    try:
        handle = open(args.csv, encoding="utf-8-sig", errors="replace", newline="")
    except OSError as exc:
        print(f"[错误] 无法读取：{exc}")
        return 2

    with handle:
        reader = csv.DictReader(handle)
        for name in (args.seq_col, args.ref_col):
            if name not in (reader.fieldnames or []):
                print(f"[错误] 列不存在：{name}。可用列：{(reader.fieldnames or [])[:20]}……")
                return 2
        for row in reader:
            seq = row.get(args.seq_col) or ""
            ref = (row.get(args.ref_col) or "").upper()
            if not seq or not ref:
                continue
            for m in PLACEHOLDER_RE.finditer(seq):
                left = rna_to_dna(seq[max(0, m.start() - args.flank):m.start()]).upper()
                right = rna_to_dna(seq[m.end():m.end() + args.flank]).upper()
                if ("-" in left or "-" in right
                        or len(left) < args.flank or len(right) < args.flank):
                    continue
                li = ref.find(left)
                if li < 0:
                    continue
                ri = ref.find(right, li + len(left))
                if ri < 0:
                    continue
                span = ri - (li + len(left))
                k = int(m.group(1))
                width = len(m.group(0))
                tested += 1
                match_k += int(span == k)
                match_width += int(span == width)
                if sampled < args.show:
                    label = (row.get("Name") or row.get("name")
                             or row.get("id") or row.get("ID") or "?")
                    samples.append((label, k, width, span))
                    sampled += 1

    if not tested:
        print("[错误] 没有任何占位符能完成锚定（检查列名与 --flank）。")
        return 2

    print("=" * 78)
    print("占位符 *[k]* 语义实证")
    print("=" * 78)
    print(f"序列列 {args.seq_col} / 参照列 {args.ref_col} / 锚点 {args.flank}bp")
    print(f"可有效锚定的占位符数: {tested:,}")
    print()
    print(f"{'记录':<34}{'k':>5}{'标记宽度':>9}{'真实跨度':>10}")
    print("-" * 78)
    for name, k, width, span in samples:
        print(f"{str(name)[:34]:<34}{k:>5}{width:>9}{span:>10}")
    print("-" * 78)
    print(f"真实跨度 == k          : {match_k:,}/{tested:,} = {100*match_k/tested:.2f}%"
          "   ← 若接近 100%，k 是「计数」")
    print(f"真实跨度 == 标记宽度   : {match_width:,}/{tested:,} = {100*match_width/tested:.2f}%"
          "   ← 若接近 100%，k 是「编号」")
    print()
    if match_k / tested > 0.95:
        print("结论：k = 该处丢失的碱基个数 → 应还原为 k 个 N。")
    elif match_width / tested > 0.95:
        print("结论：k 不是计数（真实跨度等于标记宽度）→ 不应按 k 展开。")
    else:
        print("结论：两种模型都不吻合，需人工检查该数据集的约定。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
