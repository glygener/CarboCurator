

import json
import os
import sys
from collections import defaultdict
from statistics import median

import matplotlib.pyplot as plt

# --- configuration -----------------------------------------------------------
CATS = {"disease", "gene", "drug", "species",
    "mutation", "cell_line", "cell_type",
    "DNA", "RNA"}
    
CUTOFF = {
    "disease": 0.9,
    "gene":     0.9,
    "drug":     2,
    "species": 0.9,
    "mutation": 2,
    "cell_line": 0.5,
    "cell_type": 0.5,
    "DNA": 2,
    "RNA": 2,
}

IN_DIR = "/home/cyruschauyeung/projects/CarboCurator/data/processed/bern2_ner_labels.jsonl"
OUT_DIR = "/home/cyruschauyeung/projects/CarboCurator/data/processed/bern2_ner_labels_qced.jsonl"
OUT_DIR_PLOT = "/home/cyruschauyeung/projects/CarboCurator/data/processed/ner_re/qc_plot"
BINS = 50
# -----------------------------------------------------------------------------

def bern2_qc():
    print(f"[Script] Running {__file__.split('/')[-1]}")
    probs = defaultdict(list)

    with open(IN_DIR, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            record = json.loads(line)
            for ent in record.get("Entities", []):
                obj = ent.get("obj", "")
                if obj in CATS:
            
                    val = ent.get("prob", None)
                    if isinstance(val, (int, float)) and str(val) != "nan":
                        probs[obj].append(val)

    os.makedirs(OUT_DIR_PLOT, exist_ok=True)
    print(f"[Message] Found {len(probs)} objects with probability values.")
    print("[Message] Calculating probability distribution of named entities:")

    for obj, values in probs.items():
        plt.figure()
        plt.hist(values, bins=BINS, color="teal", edgecolor="none")
        plt.xlabel("prob")
        plt.ylabel("count")
        plt.title(f"{obj} entities (n = {len(values)})")
        plt.tight_layout()
        out_file = os.path.join(OUT_DIR_PLOT, f"prob_distribution_{obj}.png")
        plt.savefig(out_file, dpi=300)
        plt.close()
        string1 = '{: <10}'.format(obj)
        string2 = '{: <11}'.format(f"N = {len(values)}")
        print(f"          {string1} {string2}", end = "", flush=True)
        print(f" (min: {min(values):.3f}, median: {median(values):.3f}, max: {max(values):.3f})")

    print("[Message] Probability distribution plots saved to:", OUT_DIR_PLOT.split("/")[-1])

def filter_entities():
    print(f"[Message] Number of lines in {IN_DIR.split('/')[-1]}: {sum(1 for _ in open(IN_DIR, 'r', encoding='utf-8'))}")
    with open(IN_DIR, "r", encoding="utf-8") as fin, \
    open(OUT_DIR, "w", encoding="utf-8") as fout:

        n_lines = 0

        for line in fin:
            if not line.strip():
                continue

            rec = json.loads(line)
            kept = []

            for ent in rec.get("Entities", []):
                obj = ent.get("obj", "")

                thresh = CUTOFF.get(obj, 0.0)
                if isinstance(ent.get("prob"), (int, float)) and ent["prob"] > thresh:
                    kept.append(ent)

            rec["Entities"] = kept
            if rec["Entities"]:
                fout.write(json.dumps(rec, ensure_ascii=False) + "\n")
            else:
                n_lines += 1
        
        print(f"[Message] Number of lines after filtering: {n_lines}")

    print(f"[Message] Filtered file saved to {OUT_DIR.split('/')[-1]}")


if __name__ == "__main__":
    bern2_qc()
    filter_entities()