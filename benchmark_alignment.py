"""量化"注释列寄存器保持"到底值多少：同一份数据，开/关 --align-col 的对比。

背景：当某个变换会改变序列长度（例如占位符 *[k]* 还原成 k 个 N），
若只改序列、不改注释列，序列与注释就会整体错位。这种错位不会报错，
只会让下游训出坏模型——所以必须量化它到底会发生多少。

做法：
    1. 用 --align-col 跑一次（注释列同步变换）
    2. 不用 --align-col 跑一次（注释列保持原样）
    3. 各自统计"序列与注释列长度不一致"的行数

用法：
    python benchmark_alignment.py <raw.csv> \
        --seq-col Sequence --align-col Structure

退出码：0 = 完成（无论结果好坏，本脚本只测量不建议）。
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))


def run_curate(csv_path, out_dir, seq_col, align_col):
    cmd = [sys.executable, os.path.join(HERE, "curate.py"), csv_path, "--out", out_dir]
    if seq_col:
        cmd += ["--seq-col", seq_col]
    if align_col:
        cmd += ["--align-col", align_col]
    result = subprocess.run(cmd, capture_output=True, text=True,
                            encoding="utf-8", errors="replace")
    if result.returncode != 0:
        raise RuntimeError(f"curate.py 失败：\n{result.stdout}\n{result.stderr}")
    cleaned = os.path.join(out_dir, "cleaned.csv")
    if not os.path.isfile(cleaned):
        raise RuntimeError("curate.py 未产出 cleaned.csv")
    return cleaned


def count_misaligned(cleaned_path, seq_col, align_col):
    """统计序列列与注释列长度不一致的行数。"""
    seq_field = f"{seq_col}_curated"
    aligned_field = f"{align_col}_curated"
    rows = misaligned = 0
    with open(cleaned_path, encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        fields = reader.fieldnames or []
        if seq_field not in fields:
            raise RuntimeError(f"产物里找不到 {seq_field}；它真的是 curate.py 的输出吗？")
        use_curated = aligned_field in fields
        if not use_curated and align_col not in fields:
            raise RuntimeError(f"产物里既没有 {aligned_field} 也没有 {align_col}")
        for row in reader:
            seq = row.get(seq_field) or ""
            other = row.get(aligned_field if use_curated else align_col) or ""
            if not seq or not other:
                continue
            rows += 1
            if len(seq) != len(other):
                misaligned += 1
    return rows, misaligned, use_curated


def main(argv=None):
    parser = argparse.ArgumentParser(description="寄存器保持的 A/B 量化")
    parser.add_argument("csv", help="原始 CSV")
    parser.add_argument("--seq-col", default=None,
                        help="序列列；默认让 curate.py 自动检测全部序列列")
    parser.add_argument("--align-col", required=True,
                        help="需要保持逐位对齐的注释列名（如 Structure）")
    parser.add_argument("--out", default=None, help="把结果写成 JSON")
    args = parser.parse_args(argv)

    tmp = tempfile.mkdtemp(prefix="align_ab_")
    try:
        dir_on = os.path.join(tmp, "on")
        dir_off = os.path.join(tmp, "off")
        print("跑配置 A（--align-col，注释列同步变换）……")
        cleaned_on = run_curate(args.csv, dir_on, args.seq_col, args.align_col)
        print("跑配置 B（不加 --align-col，注释列保持原样）……")
        cleaned_off = run_curate(args.csv, dir_off, args.seq_col, None)

        # 序列列未指定时，直接读 curate.py 记录在摘要里的配对关系，避免猜
        seq_col = args.seq_col
        if seq_col is None:
            with open(os.path.join(dir_on, "curation_summary.json"), encoding="utf-8") as f:
                summary = json.load(f)
            pairs = summary.get("align_pairs") or []
            matched = [p["sequence_column"] for p in pairs
                       if p.get("align_column") == args.align_col]
            if not matched:
                print(f"[错误] 摘要里没有注释列 {args.align_col} 的配对记录。"
                      f"请用 --seq-col 显式指定序列列。")
                return 2
            seq_col = matched[0]

        rows_on, bad_on, _ = count_misaligned(cleaned_on, seq_col, args.align_col)
        rows_off, bad_off, used_curated = count_misaligned(cleaned_off, seq_col, args.align_col)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    print()
    print("=" * 76)
    print("寄存器保持 A/B 对比")
    print("=" * 76)
    print(f"序列列 {seq_col} / 注释列 {args.align_col}")
    print()
    print(f"{'配置':<40}{'行数':>9}{'错位行':>9}{'错位率':>10}")
    print("-" * 76)
    print(f"{'A  --align-col（注释列同步变换）':<40}{rows_on:>9,}{bad_on:>9,}"
          f"{100 * bad_on / max(1, rows_on):>9.3f}%")
    label_b = "B  不加 --align-col（注释列保持原样）"
    print(f"{label_b:<40}{rows_off:>9,}{bad_off:>9,}"
          f"{100 * bad_off / max(1, rows_off):>9.3f}%")
    print("-" * 76)
    if not used_curated:
        print("说明：配置 B 的产物没有同步后的注释列，这里回退用原始注释列比较长度。")
    print(f"同步变换消除的错位行数：{bad_off - bad_on:,}")
    if rows_off and (bad_off - bad_on) > 0:
        print(f"若不做同步，约 {100 * bad_off / rows_off:.1f}% 的行会出现序列/注释错位。")
    elif bad_off == 0:
        print("本数据上未发生长度变化，故两种配置结果相同。")

    if args.out:
        payload = {
            "sequence_column": seq_col,
            "align_column": args.align_col,
            "with_align": {"rows": rows_on, "misaligned": bad_on},
            "without_align": {"rows": rows_off, "misaligned": bad_off},
            "misaligned_avoided": bad_off - bad_on,
        }
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(payload, f, ensure_ascii=False, indent=2)
        print(f"\n结果已写入: {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
