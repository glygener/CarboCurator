import pandas as pd
import re
import os
import random

abs_path = "/home/cyruschauyeung/git/CarboCurator/data_output/pubmed_abstracts_20250325"
tbl_path = "/home/cyruschauyeung/git/CarboCurator/data_output/biomarker_table/biomarker_table_20250331_4o.tsv"

SUBSET_FILE = "subset_pmids.txt"
PROGRESS_FILE = "ner_manual_validation_results.csv"

FIELDS_TO_VALIDATE = [
    "glycan_entity",
    "glycan_entity_change",
    "glycoprotein_entity_name",
    "glycoenzyme_entity_name",
    "disease_name",
    "organism",
    "specimen_type"
]

FIELD_COLORS = {
    "glycan_entity": "salmon",
    "glycan_entity_change": "green",
    "glycoprotein_entity_name": "yellow",
    "glycoenzyme_entity_name": "magenta",
    "disease_name": "cyan",
    "organism": "orange",
    "specimen_type": "plum"
}

ANSI_FIELD_COLORS = {
    "glycan_entity": "\033[91m",             
    "glycan_entity_change": "\033[92m",        
    "glycoprotein_entity_name": "\033[93m",    
    "glycoenzyme_entity_name": "\033[95m",     
    "disease_name": "\033[96m",                
    "organism": "\033[33m",                    
    "specimen_type": "\033[97m"                
}
ANSI_RESET = "\033[0m"

def load_abstracts(folder):
    abstracts = {}
    files = [os.path.join(folder, f) for f in os.listdir(folder)
             if f.startswith("pubmed_abstract") and f.endswith(".txt")]
    files.sort(key=lambda x: int(x.split("pubmed_abstract")[-1].split(".txt")[0]))
    
    for filename in files:
        try:
            with open(filename, 'r', encoding='utf-8') as f:
                text = f.read()
                match = re.search(r'PMID:\s*(\d+)', text)
                if match:
                    pmid = match.group(1)
                    abstracts[pmid] = text
                else:
                    print(f"PMID not found in {filename}")
        except Exception as e:
            print(f"Error reading {filename}: {e}")
    return abstracts

def load_biomarker_table(tsv_file):
    return pd.read_csv(tsv_file, sep="\t")

def get_entities_for_pmid(df, pmid):
    subset = df[df["pmid"].astype(str) == str(pmid)]
    entities = {field: set() for field in FIELDS_TO_VALIDATE}
    for field in FIELDS_TO_VALIDATE:
        for val in subset[field].dropna():
            for part in str(val).split("|"):
                part = part.strip()
                if part:
                    entities[field].add(part)
    return entities

def highlight_text(text, entities, color):
    for entity in entities:
        pattern = re.escape(entity)
        replacement = f'<span style="color: {color};"><strong>{entity}</strong></span>'
        text = re.sub(pattern, replacement, text)
    return text

def create_markdown_with_highlights(abstract, entities_by_field):
    highlighted_text = abstract
    for field in FIELDS_TO_VALIDATE:
        color = FIELD_COLORS.get(field, "black")
        entities = entities_by_field.get(field, set())
        highlighted_text = highlight_text(highlighted_text, entities, color)
    md_text = f"## Abstract\n\n{highlighted_text}\n"
    return md_text

def save_temp_markdown(md_text):
    filename = "temp_abstract.md"  # Always the same file name
    with open(filename, "w", encoding="utf-8") as f:
        f.write(md_text)
    return filename

def print_entities_colored(entities_by_field):
    for field in FIELDS_TO_VALIDATE:
        ansi_color = ANSI_FIELD_COLORS.get(field, "")
        entities = entities_by_field.get(field, set())
        entities_str = ", ".join(entities) if entities else "None"
        print(f"{ansi_color}{field}:{ANSI_RESET} {entities_str}")

def save_subset_pmids(sampled_pmids):
    with open(SUBSET_FILE, "w", encoding="utf-8") as f:
        for pmid in sampled_pmids:
            f.write(f"{pmid}\n")

def load_subset_pmids():
    if os.path.exists(SUBSET_FILE):
        with open(SUBSET_FILE, "r", encoding="utf-8") as f:
            return [line.strip() for line in f if line.strip()]
    return None

def save_progress(results):
    df = pd.DataFrame(results)
    df.to_csv(PROGRESS_FILE, index=False)
    print(f"Progress saved to {PROGRESS_FILE}")

def main(abs_path, tbl_path):
    abstracts_folder = abs_path
    biomarker_tsv = tbl_path
    
    print("Loading abstracts...")
    abstracts = load_abstracts(abstracts_folder)
    print(f"Loaded {len(abstracts)} abstracts.")
    
    print("Loading biomarker table...")
    df_biomarker = load_biomarker_table(biomarker_tsv)
    
    all_pmids = list(abstracts.keys())
    sample_size = 460
    subset_pmids = load_subset_pmids()
    if subset_pmids is None:
        print("No subset file found. Creating a new random subset of 460 PMIDs.")
        if len(all_pmids) < sample_size:
            print("Warning: fewer than 460 abstracts available; using all.")
            subset_pmids = all_pmids
        else:
            random.seed(42)
            subset_pmids = random.sample(all_pmids, sample_size)
        save_subset_pmids(subset_pmids)
    else:
        print(f"Loaded subset of {len(subset_pmids)} PMIDs from {SUBSET_FILE}.")
    
    if os.path.exists(PROGRESS_FILE):
        df_progress = pd.read_csv(PROGRESS_FILE)
        processed_pmids = set(df_progress["pmid"].astype(str))
        validation_results = df_progress.to_dict(orient="records")
    else:
        validation_results = []
        processed_pmids = set()
    
    remaining_pmids = [pmid for pmid in subset_pmids if pmid not in processed_pmids]
    print(f"Number of abstracts remaining for validation: {len(remaining_pmids)}")
    
    for pmid in remaining_pmids:
        print("\n==============================")
        print(f"PMID: {pmid}")
        
        abstract_text = abstracts.get(pmid, "")
        entities_by_field = get_entities_for_pmid(df_biomarker, pmid)
        
        md_text = create_markdown_with_highlights(abstract_text, entities_by_field)
        temp_filename = save_temp_markdown(md_text)
        print(f"\nThe highlighted abstract has been written to '{temp_filename}'. Please open it in your IDE for easier reading.\n")
        
        print("Recognized entities (colored):")
        print_entities_colored(entities_by_field)
        
        biomarker_rows = df_biomarker[df_biomarker["pmid"].astype(str) == str(pmid)]
        if biomarker_rows.empty:
            print("No biomarker table rows found for this PMID.")
        else:
            print("Biomarker table rows for this PMID:")
            try:
                print(biomarker_rows.to_markdown(tablefmt="grid"))
            except Exception as e:
                # Fallback if to_markdown is not available
                print(biomarker_rows.to_string(index=False))
        
        result = {"pmid": pmid}
        for field in FIELDS_TO_VALIDATE:
            prompt = f"\nFor field '{ANSI_FIELD_COLORS.get(field, '')}{field}{ANSI_RESET}', enter validation result: "
            while True:
                user_input = input(prompt).strip().lower()
                if user_input in ["tp", "tn", "fp", "fn", "pr"]:
                    result[field] = user_input
                    break
                else:
                    print("Invalid input. Please enter one of: tp, tn, fp, fn, pr.")
        
        # Ask for an additional comment for this PMID.
        comment = input("\nEnter additional comment for this PMID (or leave blank): ").strip()
        result["comment"] = comment
        
        validation_results.append(result)
        save_progress(validation_results)
        
        cont = input("\nPress Enter to continue to the next abstract, or type 'q' to quit: ").strip().lower()
        if cont == "q":
            print("Exiting and saving progress.")
            break

    final_output = "ner_manual_validation_final.csv"
    pd.DataFrame(validation_results).to_csv(final_output, index=False)
    print(f"\nFinal validation results saved to {final_output}")

if __name__ == "__main__":
    main(abs_path, tbl_path)

