"""semantic-seq-curator-agent 的对比评测基准（benchmark）。

在同一个含非标准字符的测试集上，对比五组方法的可量化表现：
  1. 白板 Agent（无任何 SKILL/知识）
  2. 简单 Prompt（一句"把非标准字符当噪声删"）
  3. seismic cleanfa（模糊码 → N 替换；非 IUPAC 报错）
  4. SNIP（含非 20AA 字符 → 整条删除）
  5. 本 Agent（semantic_seq_curator）

指标（全部 0~1，越高越好）：
  - 信息保留率 retention：9 个"应保留(有生物学意义)"字符中被正确保留的比例
  - 噪声清除率 noise_removal：2 个"应清除(测序错误)"字符中被正确掩码/移除的比例
  - 可解释率 explainability：11 个字符中带可审计理由的比例
  - 端到端正确率 accuracy：(正确保留 + 正确清除) / 11
"""

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from semantic_seq_curator.agent import SemanticSeqCuratorAgent

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# 测试集 + 金标准（should_keep=True 表示该字符有生物学意义、应保留并标注）
# ---------------------------------------------------------------------------
TESTS = [
    {"seq_id": "seq_dna", "sequence": "ACGGATNAGCCTTTTTTNTTTTTACGRX",
     "seq_type": "dna", "quality": {7: 38, 18: 8, 27: 35, 28: 30}},
    {"seq_id": "seq_prot", "sequence": "MKTAYGASGAGUGAGSGAOALBX*",
     "seq_type": "protein", "quality": {}},
    {"seq_id": "seq_artifact", "sequence": "ACGT!ACGT?ACGT",
     "seq_type": "dna", "quality": {5: 10}},
]

# (seq_id, position) -> should_keep（True=应保留, False=应清除）
GOLD = {
    ("seq_dna", 7): True,      # N 在保守基序内, 高质量 → iupac_ambiguity
    ("seq_dna", 18): False,    # N 在同聚物内, 低质量 → sequencing_error
    ("seq_dna", 27): True,     # R = purine 模糊码 → iupac_ambiguity
    ("seq_dna", 28): True,     # X = 未知碱基 → iupac_ambiguity
    ("seq_prot", 12): True,    # U = 硒代半胱氨酸 → functional_modification
    ("seq_prot", 19): True,    # O = 吡咯赖氨酸 → functional_modification
    ("seq_prot", 22): True,    # B = Asx(D/N) → iupac_ambiguity
    ("seq_prot", 23): True,    # X = 未知残基 → iupac_ambiguity
    ("seq_prot", 24): True,    # * = 终止子 → functional_modification
    ("seq_artifact", 5): False,   # ! 低质量乱码 → sequencing_error
    ("seq_artifact", 10): True,   # ? 无质量证据 → unknown_artifact(保留并人工复核)
}
TOTAL = len(GOLD)
N_KEEP = sum(1 for v in GOLD.values() if v)
N_REMOVE = sum(1 for v in GOLD.values() if not v)


# ---------------------------------------------------------------------------
# 1. 本 Agent
# ---------------------------------------------------------------------------
def run_agent():
    dec = {}
    for t in TESTS:
        a = SemanticSeqCuratorAgent("nucleic_lm_pretraining", t["seq_type"], t["quality"])
        for r in a.curate(t["seq_id"], t["sequence"]):
            keep = r.action in ("keep_annotate_uncertainty", "keep_annotate_modification",
                                "keep_annotate_gap", "flag_human_review")
            dec[(r.seq_id, r.position)] = {
                "keep": keep, "remove": not keep,
                "reason": r.reason, "has_reason": bool(r.reason.strip()),
            }
    return dec


# ---------------------------------------------------------------------------
# 2. seismic cleanfa（核酸：模糊码→N；非 IUPAC 报错/移除；蛋白不支持）
# ---------------------------------------------------------------------------
def run_cleanfa():
    AMBIG = set("RYSWKMBDHVNX")
    dec = {}
    for t in TESTS:
        if t["seq_type"] == "protein":
            # cleanfa 是核酸工具，对蛋白序列无法处理（报错），保留率计 0
            for pos, keep in GOLD.items():
                if pos[0] == t["seq_id"]:
                    dec[pos] = {"keep": False, "remove": True, "reason": "", "has_reason": False}
            continue
        seq = t["sequence"].upper()
        for i, ch in enumerate(seq):
            pos = i + 1
            if (t["seq_id"], pos) not in GOLD:
                continue
            if ch in AMBIG:
                # 模糊码统一替换为 N：若原字符就是 N 则"保留位置"，否则视为"替换(丢语义)"
                if ch == "N":
                    dec[(t["seq_id"], pos)] = {"keep": True, "remove": False,
                                               "reason": "", "has_reason": False}
                else:
                    dec[(t["seq_id"], pos)] = {"keep": False, "remove": True,
                                               "reason": "", "has_reason": False}
            else:
                # 非 IUPAC 字符：报错要求手动删除
                dec[(t["seq_id"], pos)] = {"keep": False, "remove": True,
                                           "reason": "", "has_reason": False}
    return dec


# ---------------------------------------------------------------------------
# 3. SNIP（白名单=20 标准氨基酸；含非白名单字符 → 整条删除）
# ---------------------------------------------------------------------------
def run_snip():
    AA = set("ARNDCQEGHILKMFPSTWYV")
    dec = {}
    for t in TESTS:
        seq = t["sequence"].upper()
        removed_whole = any(ch not in AA and ch not in "-." for ch in seq)
        for pos, keep in GOLD.items():
            if pos[0] != t["seq_id"]:
                continue
            dec[pos] = {"keep": False, "remove": True,  # 整条删除 → 该字符被移除
                        "reason": "", "has_reason": False}
    return dec


# ---------------------------------------------------------------------------
# 4/5. 白板 Agent / 简单 Prompt：从 subagent 返回的 JSON 加载
# ---------------------------------------------------------------------------
def load_json_method(path):
    with open(path, encoding="utf-8") as f:
        rows = json.load(f)
    dec = {}
    for row in rows:
        key = (row["seq_id"], int(row["position"]))
        action = (row.get("action") or "").lower()
        keep = "keep" in action or "保留" in action
        dec[key] = {"keep": keep, "remove": not keep,
                    "reason": row.get("reason", ""),
                    "has_reason": bool((row.get("reason") or "").strip())}
    return dec


# ---------------------------------------------------------------------------
# 指标计算
# ---------------------------------------------------------------------------
def metrics(dec):
    keep_ok = remove_ok = 0
    for key, should_keep in GOLD.items():
        d = dec.get(key)
        if d is None:
            continue  # 未处理，视为错误
        if should_keep and d["keep"]:
            keep_ok += 1
        if (not should_keep) and d["remove"]:
            remove_ok += 1
    explain = sum(1 for key in GOLD if dec.get(key) and dec[key]["has_reason"])
    return {
        "retention": keep_ok / N_KEEP,
        "noise_removal": remove_ok / N_REMOVE,
        "explainability": explain / TOTAL,
        "accuracy": (keep_ok + remove_ok) / TOTAL,
    }


def main():
    methods = {
        "白板 Agent": load_json_method(os.path.join(HERE, "benchmark", "blank_agent.json")),
        "简单 Prompt": load_json_method(os.path.join(HERE, "benchmark", "simple_prompt.json")),
        "seismic cleanfa": run_cleanfa(),
        "SNIP": run_snip(),
        "本 Agent": run_agent(),
    }
    results = {name: metrics(dec) for name, dec in methods.items()}

    out_dir = os.path.join(HERE, "benchmark")
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "results.json"), "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=2)

    print("=" * 78)
    print("对比评测结果（金标准 11 个非标准字符：应保留 9 / 应清除 2）")
    print("=" * 78)
    print(f"{'方法':<16}{'信息保留率':>12}{'噪声清除率':>12}{'可解释率':>12}{'端到端正确率':>12}")
    print("-" * 64)
    for name, m in results.items():
        print(f"{name:<16}{m['retention']*100:>11.0f}%{m['noise_removal']*100:>11.0f}%"
              f"{m['explainability']*100:>11.0f}%{m['accuracy']*100:>11.0f}%")
    print("=" * 78)
    return results


if __name__ == "__main__":
    main()
