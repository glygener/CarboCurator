# Script 1: preprocess_get_pmid.py
# Description: This script fetches PubMed IDs (PMIDs) based on a search query using the NCBI E-utilities API.
# Usage: python proprocess_get_pmid.py
# Input: A search query entered through the script below.
# Output: A list of PMIDs saved to a text file.
# Dependencies: Biopython
# Note: Set your NCBI API key in the environment variable NCBI_API_KEY for higher rate limits.
# WARNING: If the number of records in a year exceeds 9999, the script will fail to retrieve all records.

__version__ = "1.0.0"


from Bio import Entrez
from datetime import datetime
import time
import os

# User-defined parameters #######################################################
term = "(((glycomics[Title/Abstract]) OR (polysaccharides[MeSH Terms])) AND ((biomarkers[MeSH Terms])) OR (disease[MeSH Terms]))"
start_year = 1900
end_year = "current" # Use "current" to get the current year
Entrez.email = "your.email@here.com" # Optional
NCBI_API_KEY = "YOUR_API_KEY" # Optional: replace with your key or set env var
output_translation_dir = "/home/cyruschauyeung/projects/CarboCurator/data/processed/corpus/translaton_NER.txt"
output_pmid_list_dir = "/home/cyruschauyeung/projects/CarboCurator/data/processed/corpus/pmid_NER.txt"
#################################################################################

# System-defined parameters #####################################################
API_LIMIT = 9999 # retstart won't work if records > 9999. Trust me
if end_year == "current":
    end_year = datetime.now().year
if start_year > end_year:
    raise ValueError(f"Start year {start_year} is greater than end year {end_year}.")
#################################################################################

def check_api_key():
    NCBI_API_KEY = os.getenv("NCBI_API_KEY")
    if NCBI_API_KEY != "YOUR_API_KEY":
        Entrez.api_key = NCBI_API_KEY
        print("[Message] Environmental variable NCBI_API_KEY detected and using.")
    else:
        print("[Message] No NCBI_API_KEY detected. Using default rate limits.")
        print("          Consider setting NCBI_API_KEY for higher rate limits.")
        print("          See https://support.nlm.nih.gov/kbArticle/?pn=KA-05317 for more information.")

def fetch_translation(term, start_year, end_year):
    esearch_params = {"db":"pubmed", "term":term,
                      "mindate":f"{start_year}/01/01",
                      "maxdate":f"{end_year}/12/31",
                      "retmax":0}
    handle = Entrez.esearch(**esearch_params)
    translation = Entrez.read(handle)['QueryTranslation']
    handle = Entrez.esearch(**esearch_params) # handle closes after reading, that's why
    count = Entrez.read(handle)['Count']
    handle.close()
    return translation, count

def esearch_ids(term, a, b, retmax=API_LIMIT):
    h = Entrez.esearch(
        db="pubmed", term=term,
        mindate=f"{a}/01/01",
        maxdate=f"{b}/12/31",
        usehistory="n",
        retmax=retmax,
        sort="pub_date")
    r = Entrez.read(h)
    h.close()
    return int(r["Count"]), r["IdList"]

def batch_retrieval(term, a, b): # recursive splitting, can be improved
    queue = [(a, b)]
    pmids = []

    while queue:
        a, b = queue.pop()
        count, ids = esearch_ids(term, a, b)
        print(f"{a}-{b}(N={count})", end = " ", flush = True)

        if count == 0:
            continue
        if count <= API_LIMIT or a == b: # leaf
            pmids.extend(ids)
            if count > len(ids):
                print(f"\n[Warning] {a}-{b} clipped at {API_LIMIT}", file=sys.stderr)
        else:
            mid = (a + b) // 2
            queue.append((a, mid))
            queue.append((mid + 1, b))
        time.sleep(0.4)

    return pmids

def write_translation(translation):
    with open(output_translation_dir, "w") as f:
        f.write(translation)
        f.write("\n")


def write_pmid_list(pmids):
    with open(output_pmid_list_dir, "w") as f:
        f.write("\n".join(pmids))
        f.write("\n")

def main():
    print(f"[Script] Running {__file__.split('/')[-1]}")
    check_api_key()
    print(f"[Message] Search term: {term}")
    esearch_args = (term, start_year, end_year)

    # Get Query Translation
    print("[Message] (1) Getting the query translation of the search term...")
    translation, count = fetch_translation(*esearch_args)
    write_translation(translation)
    print(f"[Message] Saved translation to {output_translation_dir.split('/')[-1]}")
    print(f"[Message] Total records found: {count}")

    # Get PMID List
    print("[Message] (2) Fetching PMIDs...")
    print("[Message] Recursive splitting year range: ", end = " ", flush = True)
    pmids = batch_retrieval(*esearch_args)
    write_pmid_list(pmids)
    print(f"\n[Message] Saved {len(pmids)} PMIDs to {output_pmid_list_dir.split('/')[-1]}")
    print("[Message] Done.")


if __name__ == "__main__":
    main()
