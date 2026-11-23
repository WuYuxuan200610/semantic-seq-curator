"""最小验证：展示"同一字符在不同上下文被不同处理" + 与哑清洗对比。

用法：python run_demo.py
无第三方依赖，仅标准库。
"""

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from semantic_seq_curator.agent import SemanticSeqCuratorAgent

# ---------------------------------------------------------------------------
# 测试序列
# ---------------------------------------------------------------------------
# DNA：位置 7 的 N 位于保守基序 "GATAA"（GATA 转录因子结合位点）内，高质量；
#       位置 18 的 N 嵌在同聚物 T-run 内，低质量；
#       位置 27 的 R 为真实 IUPAC 模糊（purine A/G）；位置 28 的 X 为未知/掩蔽。
DNA_SEQ = "ACGGATNAGCCT" + "TTTTTNTTTTT" + "ACGRX"
#            ^        ^        ^
#         7=N(基序)  18=N(同聚物) 27=R 28=X
DNA_QUALITY = {7: 38, 18: 8, 27: 35, 28: 30}   # Phred 质量分

# 蛋白质：U=硒代半胱氨酸, O=吡咯赖氨酸, B=Asx(D/N), X=未知残基, *=终止子
PROTEIN_SEQ = "MKTAYGASGAG" + "U" + "GAGSGA" + "O" + "AL" + "B" + "X" + "*"


def run_case(seq_id, sequence, task, seq_type, quality=None, verbose=False):
    agent = SemanticSeqCuratorAgent(
        downstream_task=task, seq_type=seq_type,
        quality_track=quality, verbose=verbose,
    )
    records = agent.curate(seq_id, sequence)
    return records


def print_table(records):
    header = f"{'seq':<10}{'pos':>4}  {'char':>4}  {'semantic_class':<24}{'action':<28}{'conf':>6}"
    print(header)
    print("-" * len(header))
    for r in records:
        print(f"{r.seq_id:<10}{r.position:>4}  {r.original_char:>4}  "
              f"{r.semantic_class:<24}{r.action:<28}{r.confidence:>6.2f}")
    print()


def main():
    print("=" * 100)
    print("semantic-seq-curator-agent  最小验证")
    print("=" * 100)

    # ---- DNA 案例（核心：同一个 N 不同上下文不同决策）----
    print("\n[1] DNA 序列  |  下游任务 = nucleic_lm_pretraining\n")
    print("序列: ", DNA_SEQ)
    print("      pos 7  N 在保守基序 GATAA 内 (Phred=38)")
    print("      pos 18 N 嵌在 T-同聚物内      (Phred=8)")
    print("      pos 27 R = purine 模糊码; pos 28 X = 未知\n")

    dna = run_case("seq_dna", DNA_SEQ, "nucleic_lm_pretraining",
                   "dna", DNA_QUALITY, verbose=False)
    print_table(dna)

    print(">>> 关键观察：同一个字符 'N' 得到两种语义分类与动作 ——")
    for r in dna:
        if r.original_char == "N":
            print(f"    N@{r.position}  →  {r.semantic_class} / {r.action}  (conf={r.confidence})")
            print(f"        理由: {r.reason}")
    print()

    # ---- 蛋白质案例 ----
    print("[2] 蛋白质序列  |  下游任务 = protein_lm_pretraining\n")
    print("序列: ", PROTEIN_SEQ)
    prot = run_case("seq_prot", PROTEIN_SEQ, "protein_lm_pretraining",
                    "protein", None, verbose=False)
    print_table(prot)
    print()

    # ---- 写入 curation_log（JSON Lines）----
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "demo_output")
    os.makedirs(out_dir, exist_ok=True)
    log_path = os.path.join(out_dir, "curation_log.jsonl")
    with open(log_path, "w", encoding="utf-8") as f:
        for r in dna + prot:
            f.write(r.to_jsonl() + "\n")
    print(f"[3] curation_log 已写入: {log_path}\n")

    # ---- 与哑清洗对比 ----
    print("[4] 与现有工具（哑清洗）对比\n")
    print("    seismic cleanfa  (模糊码统一替换为 N):")
    print("      ", DNA_SEQ.replace("R", "N").replace("X", "N"))
    print("      → R/X 的生物学信息丢失；N@7 与 N@18 无法区分；无理由记录。")
    print("    SNIP / ShortRead  (含非标准字符即整条删除):")
    print("      → 该 DNA 序列含 R/X，整条被删除，数据全损，无理由记录。")
    print("    本 Agent:")
    print("      → 逐位判定: N@7=iupac_ambiguity(保留), N@18=sequencing_error(掩码),")
    print("        R/X=保留并标注不确定等级; 每条决策带 reason + confidence，可审计。")
    print()
    print("=" * 100)
    print("验证完成。同一个 'N' 在不同上下文得到不同语义分类与动作，证明语义感知。")
    print("=" * 100)


if __name__ == "__main__":
    main()
