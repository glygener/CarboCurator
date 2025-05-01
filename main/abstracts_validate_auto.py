import pandas as pd
import re
import os

def load_abstracts(folder):
    """
    Load abstracts from all text files in the given folder and build a dictionary keyed by PMID.
    Assumes files are named 'pubmed_abstract1.txt', 'pubmed_abstract2.txt', ..., with no leading zeros.
    """
    abstracts = {}
    # List all files in the folder that match the pattern.
    files = [os.path.join(folder, f) for f in os.listdir(folder) if f.startswith("pubmed_abstract") and f.endswith(".txt")]
    # Sort files numerically by the number following 'pubmed_abstract'
    files.sort(key=lambda x: int(x.split("pubmed_abstract")[-1].split(".txt")[0]))
    
    for filename in files:
        with open(filename, 'r', encoding='utf-8') as f:
            text = f.read()
            # Extract the PMID from the abstract text using regex.
            match = re.search(r'PMID:\s*(\d+)', text)
            if match:
                pmid = match.group(1)
                abstracts[pmid] = text
            else:
                print(f"PMID not found in {filename}")
    return abstracts

def validate_sentence(abstracts, sentence, pmid):
    """
    Return True if the extracted sentence is found in the corresponding abstract for the given PMID.
    """
    if pd.isnull(sentence):
        return False
    abstract = abstracts.get(str(pmid), "")
    # Ensure sentence is a string.
    sentence = str(sentence)
    return sentence in abstract

def main():
    # Define folder containing all abstract files and the TSV file with NER output.
    folder = "/home/cyruschauyeung/git/CarboCurator/data_output/pubmed_abstracts_20250325"
    tsv_file = "/home/cyruschauyeung/git/CarboCurator/data_output/biomarker_table/biomarker_table_20250331_4o.tsv"
    
    # Load all abstracts into a dictionary keyed by PMID.
    abstracts = load_abstracts(folder)
    
    # Read the NER output TSV file into a pandas DataFrame.
    df = pd.read_csv(tsv_file, sep="\t")
    
    # Validate each row by checking whether the extracted sentence is present in the corresponding abstract.
    df["valid"] = df.apply(lambda row: validate_sentence(abstracts, row["sentence_in_abstract"], row["pmid"]), axis=1)
    
    # Report any rows that did not validate.
    invalid_rows = df[~df["valid"]]
    if not invalid_rows.empty:
        print("The following rows did not match their abstracts:")
        print(invalid_rows[["pmid", "sentence_in_abstract"]])
    else:
        print("All rows validated successfully!")
    
    # Save the full validation results to a CSV file for further review.
    output_file = "validation_results_all.csv"
    df.to_csv(output_file, index=False)
    print(f"Validation results saved to {output_file}")

if __name__ == "__main__":
    main()
