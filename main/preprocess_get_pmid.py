# Script 1: pmid_fetch.py
# Description: This script fetches PubMed IDs (PMIDs) based on a search query using the NCBI E-utilities API.
# Usage: python pmid_fetch.py
# Input: A search query entered through the script below.
# Output: A list of PMIDs saved to a text file.
# Dependencies: Biopython
# Note: Set your NCBI API key in the environment variable NCBI_API_KEY for higher rate limits.

from Bio import Entrez
from datetime import datetime
import time
import os

# User-defined parameters #######################################################
term = "Glycan biomarker",
start_year = 1900
end_year = "current" # Use "current" to get the current year
Entrez.email = "your.email@here.com" # Replace with your email address
NCBI_API_KEY = "YOUR_API_KEY" # Optional: replace with your key or set env var
#################################################################################

def check_api_key():
    NCBI_API_KEY = os.getenv("NCBI_API_KEY")
    if NCBI_API_KEY is not "YOUR_API_KEY":
        Entrez.api_key = NCBI_API_KEY
        print("[Message] Environmental variable NCBI_API_KEY detected and using.")
    else:
        print("[Message] No NCBI_API_KEY detected. Using default rate limits.")
        print("          Consider setting NCBI_API_KEY for higher rate limits.")
        print("          See https://support.nlm.nih.gov/kbArticle/?pn=KA-05317 for more information.")

def count_pmids(term, start_year, end_year):
    print(f"{start_year}-{end_year}", end = " ", flush = True)
    handle = Entrez.esearch(
        db="pubmed", term=term,
        mindate=f"{start_year}/01/01",
        maxdate=f"{end_year}/12/31",
        retmax=0)
    count = int(Entrez.read(handle)["Count"])
    if (end_year - start_year == 1):
        if count > API_LIMIT:
            print(f"[Warning] There are {API_LIMIT} hits in {start_year} - {end_year}. Not all will be retrieved.")
    return count
    handle.close()

def fetch_translation(term, start_year, end_year):
    handle = Entrez.esearch(
        db="pubmed", term=term,
        mindate=f"{start_year}/01/01",
        maxdate=f"{end_year}/12/31",
        retmax=0)
    translation = Entrez.read(handle)['QueryTranslation']
    return translation
    handle.close()

def fetch_pmids(term, start_year, end_year, batch=API_LIMIT):
    handle = Entrez.esearch(
        db="pubmed", term=term,
        mindate=f"{start_year}/01/01",
        maxdate=f"{end_year}/12/31",
        retmax=batch,
        sort="pub_date")
    return Entrez.read(handle)["IdList"]
    handle.close()


def batch_retrieval_by_recursive_splitting(term, start_year, end_year):
    queue = [(start_year, end_year)]
    all_pmids = []

    while queue:
        #print(f"Queue: {queue}")
        a, b = queue.pop()
        n = count_pmids(term, a, b)
        if n == 0:
            continue # nothing in this window
        if n <= API_LIMIT or a == b: # leaf – safe to fetch
            all_pmids.extend(fetch_pmids(term, a, b))
            time.sleep(0.3)
        else: # split the busier half
            mid = (a + b) // 2
            queue.append((a, mid))
            queue.append((mid + 1, b))
    return all_pmids

##################################################################################

if __name__ == "__main__":
    # System variables
    API_LIMIT = 9999 # All because of this...
    if end_year == "current":
        end_year = datetime.now().year
    common_args = (term, start_year, end_year)

    # Get Query Translation
    print("[Message] (1) Getting the query translation of the search term...")
    translation = fetch_translation(*common_args)
    with open("translation.txt", "w") as f:
        f.write(translation)
        f.write("\n")
    print(f"[Message] Saved translation to translation.txt")

    # Get PMID List
    print("[Message] (2) Fetching PMIDs...")
    print("[Message] Recursive splitting year range: ", end = " ", flush = True)
    pmids = batch_retrieval_by_recursive_splitting(*common_args)
    with open("pmid_list.txt", "w") as f:
        f.write("\n".join(pmids))
        f.write("\n")
    print(f"\n[Message] Saved {len(pmids):,} PMIDs to pmid_list.txt")
    print("[Message] Done.")
