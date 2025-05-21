# Script 2: preprocess_get_abst.py
# Description: This script fetches abstracts from PubMed using the Entrez API.
# Usage: python preprocess_get_abst.py
# Input: A list of PubMed IDs separated by new lines, in the data_input directory.
# Output: A TSV file and text files containing the fetched abstracts, in the data_output directory.
# Dependencies: Biopython, pandas, tqdm
# Note: Adjust default batch size if needed.

import os
import pandas as pd
from tqdm import tqdm
from Bio import Entrez

# User‑defined parameters ######################################################
INPUT_FILE   = "pmid_list.txt"              # PMID list located in data_input/
OUTPUT_TAG   = "abstracts"                  # Base name for output directory / TSV
BATCH_SIZE   = 50                           # PMIDs fetched per batch (≤ 1000)
EMAIL        = "your.username@example.com"  # Contact email for NCBI
NCBI_API_KEY = "YOUR_API_KEY"               # Or set the env var NCBI_API_KEY
################################################################################

def check_api_key():
    api_key_env = os.getenv("NCBI_API_KEY")
    if api_key_env and api_key_env != "YOUR_API_KEY":
        Entrez.api_key = api_key_env
        print("[Message] Environmental variable NCBI_API_KEY detected and using.")
    elif NCBI_API_KEY and NCBI_API_KEY != "YOUR_API_KEY":
        Entrez.api_key = NCBI_API_KEY
        print("[Message] Using NCBI_API_KEY from user parameters.")
    else:
        print("[Message] No NCBI_API_KEY detected. Using default rate limits.")
        print("          Consider setting NCBI_API_KEY for higher rate limits.")
        print("          See https://support.nlm.nih.gov/kbArticle/?pn=KA-05317 for more information.")


def fetch_xml(pmids):
    handle = Entrez.efetch(db="pubmed", id=",".join(pmids), retmode="xml")
    xml = Entrez.read(handle)
    handle.close()
    return xml


def extract_fields(article):
    medline_citation = article["MedlineCitation"]
    pmid = medline_citation["PMID"]

    article_data = medline_citation.get("Article", {})
    title = article_data.get("ArticleTitle", "")
    abstract = " ".join(article_data.get("Abstract", {}).get("AbstractText", ""))
    
    # Extract keywords
    keyword_list = ""
    if "KeywordList" in medline_citation and medline_citation["KeywordList"]:
        keywords_field = medline_citation["KeywordList"]
        if isinstance(keywords_field[0], list):
            keywords = keywords_field[0]
        else:
            keywords = keywords_field
        keyword_list = "|".join(keywords) if keywords else ""

    # Extract chemical substances
    chemical_names, chemical_rn, chemical_uid = [], [], []
    
    chemical_list = medline_citation.get("ChemicalList", [])
    for chem in chemical_list:
        chemical_names.append(chem["NameOfSubstance"])
        chemical_rn.append(chem.get("RegistryNumber", ""))
        chemical_uid.append(chem["NameOfSubstance"].attributes.get("UI", ""))
    
    chemical_names = "|".join(chemical_names) if chemical_names else ""
    chemical_rn = "|".join(chemical_rn) if chemical_rn else ""
    chemical_uid = "|".join(chemical_uid) if chemical_uid else ""
    
    # Extract MeSH headings
    mesh_names, mesh_uid = [], []
    
    mesh_list = medline_citation.get("MeshHeadingList", [])
    for mesh in mesh_list:
        mesh_names.append(mesh["DescriptorName"])
        mesh_uid.append(mesh["DescriptorName"].attributes.get("UI", ""))
    
    mesh_names = "|".join(mesh_names) if mesh_names else ""
    mesh_uid = "|".join(mesh_uid) if mesh_uid else ""
    
    return {
        "PMID": int(pmid),
        "Title": title,
        "Abstract": abstract,
        "Keywords": keyword_list,
        "Chemical_names": chemical_names,
        "Chemical_rn": chemical_rn,
        "Chemical_uid": chemical_uid,
        "MeSH_names": mesh_names,
        "MeSH_uid": mesh_uid
    }

def save_tsv(data, output_path):
    df = pd.DataFrame(data)
    df.to_csv(output_path, sep='\t', index=False)


def save_txt(output_dir, extracted_fields, row):
    row = str(row).zfill(6)
    abs_filename = os.path.join(output_dir, f"abstract_{row}.txt") # abstracts_curate.py dependency
    with open(abs_filename, "w", encoding="utf-8") as abs_file:   # abs = abstract
        abs_file.write(
            f"PMID:\t{extracted_fields['PMID']}\n"
            f"Title:\t{extracted_fields['Title']}\n"
            f"Abstract:\t{extracted_fields['Abstract']}\n\n"
            f"Keywords:\t{extracted_fields['Keywords']}\n"
            f"Chemical_names:\t{extracted_fields['Chemical_names']}\n"
            f"Chemical_rn:\t{extracted_fields['Chemical_rn']}\n"
            f"Chemical_uid:\t{extracted_fields['Chemical_uid']}\n"
            f"MeSH_names:\t{extracted_fields['MeSH_names']}\n"
            f"MeSH_uid:\t{extracted_fields['MeSH_uid']}\n\n"
            f"Processing_id:\t{row}\n"
        )

if __name__ == "__main__":
    # Paths
    script_dir  = os.path.dirname(os.path.realpath(__file__))
    input_path  = os.path.join(script_dir, "..", "data_input", INPUT_FILE)
    output_dir  = os.path.join(script_dir, "..", "data_output", OUTPUT_TAG)  # For storing TXT abstracts
    output_path = os.path.join(script_dir, "..", "data_output", OUTPUT_TAG + ".tsv")
    os.makedirs(output_dir, exist_ok=True)

    Entrez.email = EMAIL
    check_api_key()

    # Main
    try:
        with open(input_path, 'r') as file:
            pmid_list = [pmid.strip() for pmid in file]
        print(f"[\033[92mProcess Started\033[0m] Batch size = {BATCH_SIZE}. Fetching {len(pmid_list)} PMIDs... ")

    except FileNotFoundError:
        print(f"[\033[91mError\033[0m] File not found: {input_path}")
        exit(0)

    tsv = []
    row = 1  # Index; row number in TSV data.

    for i in tqdm(range(0, len(pmid_list), BATCH_SIZE),
                  desc="Processing", unit="batch", colour="green"):
        batch_pmids = pmid_list[i:i+BATCH_SIZE]
        xml = fetch_xml(batch_pmids)

        for article in xml["PubmedArticle"]:
            extracted_fields = extract_fields(article)
            save_txt(output_dir, extracted_fields, row)
            processing_id = str(row).zfill(6)
            tsv.append({"Processing_id": processing_id} | extracted_fields)
            row += 1

    save_tsv(tsv, output_path)
    print(f"[\033[92mComplete\033[0m] Abstracts fetched and parsed to TSV and TXT.\n"
          f" Text files saved to:\t{output_dir}\n TSV table saved to:\t{output_path}")