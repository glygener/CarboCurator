# Script 3: preprocess_qc_abst.py
# Description: This script performs quality control on retrievedPubMed records.
# Usage: python preprocess_qc_abst.py
# Input: A JSONL file containing PubMed records.
# Output: A JSONL file with duplicate records removed and quality control checks performed.
# Dependencies: NA
# Note: NA
# Warning: NA

__version__ = "1.0.0"

import json

# User-defined parameters #######################################################
input_file = "/home/cyruschauyeung/projects/CarboCurator/data/processed/pubmed_records.jsonl"
output_file = "/home/cyruschauyeung/projects/CarboCurator/data/processed/pubmed_records_qced.jsonl"
#################################################################################


def deduplicate(records):
    seen = {}
    deduped = []
    dupe_index = []
    for rec in records:
        key = rec["Title"].strip()
        if key not in seen:
            seen[str(key)] = True
            deduped.append(rec)
        else:
            dupe_index.append(rec["Processing_id"])
    if len(dupe_index) > 0:
        print(f"[Warning] Found {len(dupe_index)} duplicate records:")
        print(f"[Warning] Processing IDs: {', '.join(dupe_index)}")
    else:
        print(f"[Message] No duplicate records found (after removing empty abstracts).")
            
    return deduped

def remove_lines_with_empty_abstract(records):
    filtered_records = [rec for rec in records if rec["With_abstract"] == True]
    print(f"[Message] Removed {len(records) - len(filtered_records)} records with empty abstracts.")
    return filtered_records

def quality_control(records, casefold: bool = False):

    no_abstract = sum(1 for rec in records if not rec["With_abstract"])
    totala = len(records)
    percenta = (no_abstract / totala) * 100
    print(f"[Message] {no_abstract}/{totala} articles ({percenta:.2f}%) have an empty abstract.")


    review_count = sum(1 for rec in records if rec["Is_review"])
    totalr = len(records)
    percentr = (review_count / totalr) * 100
    print(f"[Message] {review_count}/{totalr} articles ({percentr:.2f}%) contain the string 'review'.")

def main():
    print(f"[Script] Running {__file__.split('/')[-1]}")
    with open(input_file, "r", encoding="utf-8") as f:
        records = [json.loads(line) for line in f]
    print(f"[Message] Loaded {len(records)} records from {input_file.split('/')[-1]}")

    quality_control(records)

    filtered_records = remove_lines_with_empty_abstract(records)

    original_len = len(filtered_records)
    records = deduplicate(filtered_records)
    if len(records) < original_len:
        print(f"[Message] Removed {original_len - len(records)} duplicate records.")

    with open(output_file, "w", encoding="utf-8") as f:
        for rec in records:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")

    print(f"[Message] Done. Saved {len(records)} records to {output_file.split('/')[-1]}")

if __name__ == "__main__":
    main()