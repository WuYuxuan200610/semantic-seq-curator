"""独立复核 curate.py 的产出：输出里是否真的没有"喂不进模型"的字符。

故意不复用 curate.py 的任何代码——它信任工具自己的报告就失去意义了。
直接重新读一遍成品，逐字符检查。

用法：
    python validate_model_ready.py <cleaned.csv>

退出码：0 = 全部通过；1 = 发现问题。
"""
from __future__ import annotations

import argparse
import collections
import csv
import re
import sys

# 与 curate.py 中的字母表定义保持一致的独立副本（故意重复，避免共用实现）
PROFILES = {
    "dna": set("ACGTN"),
    "rna": set("ACGUN"),
    "nucleic-mixed": set("ACGTUN"),
    "dna-iupac": set("ACGTNRYWSKMBDHV"),
    "rna-iupac": set("ACGUNRYWSKMBDHV"),
    "protein": set("ACDEFGHIKLMNPQRSTVWYX"),
    "protein-extended": set("ACDEFGHIKLMNPQRSTVWYXBZJUO*"),
}
# 这些字符一旦出现在序列列里，就一定喂不进任何核酸/蛋白模型
FORBIDDEN = set("-.*[]/\\ \t\r\n") | set("abcdefghijklmnopqrstuvwxyz")
PLACEHOLDER_RE = re.compile(r"\*\[ *(\d+)\]\*")


def main(argv=None):
    parser = argparse.ArgumentParser(description="独立复核模型就绪产物")
    parser.add_argument("cleaned_csv", help="curate.py 生成的 cleaned.csv")
    parser.add_argument("--align-col", action="append", default=None,
                        help="需与序列等长、逐位对齐的注释列，格式 序列列=注释列"
                             "（如 rna_seq=structure_curated）；可重复。"
                             "不参与字母表检查，只检查长度")
    args = parser.parse_args(argv)

    rows = checked = bad_char = bad_len = bad_delta = 0
    hist = collections.Counter()
    seq_cols = []
    align_pairs = []
    for spec in (args.align_col or []):
        if "=" not in spec:
            print(f"[错误] --align-col 需要 序列列=注释列 格式，收到：{spec}")
            return 2
        seq_part, align_part = (p.strip() for p in spec.split("=", 1))
        align_pairs.append((seq_part, align_part))

    try:
        handle = open(args.cleaned_csv, encoding="utf-8-sig", newline="")
    except OSError as exc:
        print(f"[错误] 无法读取：{exc}")
        return 2

    with handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        seq_cols = [c[:-len("_curated")] for c in fieldnames
                    if c.endswith("_curated") and f"{c[:-len('_curated')]}_alphabet" in fieldnames]
        if not seq_cols:
            print("[错误] 找不到任何 '<列>_curated' + '<列>_alphabet' 配对；"
                  "这真的是 curate.py 的产出吗？")
            return 2
        for row in reader:
            rows += 1
            for col in seq_cols:
                value = row.get(f"{col}_curated") or ""
                if not value:
                    continue
                checked += len(value)
                hist.update(value)
                alphabet = row.get(f"{col}_alphabet") or ""
                allowed = PROFILES.get(alphabet, set())
                bad_char += sum(1 for ch in value if ch not in allowed)
            for seq_part, align_part in align_pairs:
                a = row.get(align_part) or ""
                b = row.get(f"{seq_part}_curated") or ""
                if a and b and len(a) != len(b):
                    bad_len += 1

    print("=" * 74)
    print("独立复核：模型就绪产物")
    print("=" * 74)
    print(f"文件                      : {args.cleaned_csv}")
    print(f"行数                      : {rows:,}")
    print(f"序列列                    : {seq_cols}")
    print(f"检查字符数                : {checked:,}")
    print(f"字母表外字符              : {bad_char}   (必须为 0)")
    print(f"出现过的禁用字符(应无)    : {sorted(set(hist) & FORBIDDEN)}")
    print(f"实际出现的字符全集        : {''.join(sorted(hist))}")
    if align_pairs:
        print(f"对齐列长度不符的行        : {bad_len}   (必须为 0)")
    ok = bad_char == 0 and bad_len == 0
    print()
    print("结论:", "全部通过 ✅ 输出可直接进模型输入口" if ok else "存在问题 ❌")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
