# Script 2: preprocess_get_abst.py
# Description: This script fetches PubMed records based on a list of PMIDs and extracts relevant information.
# Usage: python preprocess_get_abst.py
# Input: A list of PMIDs from a text file.
# Output: A JSONL file containing the extracted information.
# Dependencies: biopython, beautifulsoup4, tqdm
# Note: Set your NCBI API key in the environment variable NCBI_API_KEY for higher rate limits.
# Warning: NA

__version__ = "1.0.0"

import os, re, html, json, pathlib
from itertools import islice
from typing import Iterable, List
from collections import defaultdict
from Bio import Entrez
from bs4 import BeautifulSoup
from tqdm import tqdm

# User-defined parameters #######################################################
PMID_FILE       = "/home/cyruschauyeung/projects/CarboCurator/data_input/pmid_list.txt"          # one PMID per line
OUT_FILE        = "/home/cyruschauyeung/projects/CarboCurator/data_output/pubmed_test.jsonl"
BATCH_SIZE      = 100                  # try 100
Entrez.email    = "YOUR.EMAIL@example.com"
NCBI_API_KEY    = "YOUR_API_KEY"       # or set env var NCBI_API_KEY
#################################################################################


def check_api_key():
    env_key = os.getenv("NCBI_API_KEY")
    if env_key and env_key != "YOUR_API_KEY":
        Entrez.api_key = env_key
        print("[Message] Using NCBI_API_KEY from environment.")
    elif NCBI_API_KEY and NCBI_API_KEY != "YOUR_API_KEY":
        Entrez.api_key = NCBI_API_KEY
        print("[Message] Using NCBI_API_KEY from script settings.")
    else:
        print("[Message] No NCBI_API_KEY detected. Using default rate limits.\n"
              "          See https://support.nlm.nih.gov/kbArticle/?pn=KA-05317 "
              "for more information.")


def chunked(iterable: Iterable[str], size: int) -> Iterable[List[str]]:
    it = iter(iterable)
    while chunk := list(islice(it, size)):
        yield chunk


def fetch_xml(pmids: List[str]) -> str:
    try:
        with Entrez.efetch(db="pubmed", id=",".join(pmids), retmode="xml") as h:
            return h.read()
    except Exception as e:
        print(f"[Error] Failed to fetch records: {e}")
        print("[Message] Adjusting batch size and retrying...")
        pmids = pmids[:len(pmids)//2]
        time.sleep(5)
        return fetch_xml(pmids)
    
    with Entrez.efetch(db="pubmed", id=",".join(pmids), retmode="xml") as h:
        return h.read()


def clean(text: str) -> str:
    text = html.unescape(text)
    text = text.replace("\xa0", " ")
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def parse_pubmed_xml(xml_raw: str, start_pid: int) -> List[dict]:
    soup = BeautifulSoup(xml_raw, "lxml-xml")
    records, pid = [], start_pid

    for art in soup.find_all("PubmedArticle"):
        pmid = int(art.PMID.text.strip())
        title = clean(art.find("ArticleTitle").get_text(" ", strip=True))
        abstract = " ".join(
            clean(t.get_text(" ", strip=True)) for t in art.find_all("AbstractText")
        )
        keywords = [clean(k.get_text(strip=True)) for k in art.find_all("Keyword")]

        chem_names, chem_rn, chem_uid = [], [], []
        for chem in art.find_all("Chemical"):
            name_tag = chem.find("NameOfSubstance")
            chem_names.append(name_tag.get_text(strip=True))
            chem_rn.append(chem.RegistryNumber.get_text(strip=True))
            chem_uid.append(name_tag["UI"])

        mesh_names, mesh_uid = [], []
        for mh in art.find_all("MeshHeading"):
            desc = mh.find("DescriptorName")
            mesh_names.append(desc.get_text(strip=True))
            mesh_uid.append(desc["UI"])
        
        # Label articles without abstracts
        if abstract == "":
            With_abstract = False
        else:
            With_abstract = True

        # Label articles with the word "review"
        if re.search(r"\breview\b", title, re.I) or re.search(r"\breview\b", abstract, re.I):
            Is_review = True
        else:
            Is_review = False

        records.append({
            "Processing_id": str(pid).zfill(6),
            "With_abstract": With_abstract,
            "Is_review": Is_review,
            "PMID": pmid,
            "Title": title,
            "Abstract": abstract,
            "Keywords": keywords,
            "Chemical_names": chem_names,
            "Chemical_rn": chem_rn,
            "Chemical_uid": chem_uid,
            "MeSH_names": mesh_names,
            "MeSH_uid": mesh_uid,
        })
        pid += 1

    return records


def main():
    print(f"[Script] Running {__file__.split('/')[-1]}")
    check_api_key()

    pmids = [ln.strip() for ln in pathlib.Path(PMID_FILE).read_text().splitlines()
             if ln.strip()]
    next_pid = 1

    with open(OUT_FILE, "w", encoding="utf-8") as out, \
         tqdm(total=len(pmids), colour="green",
              desc="[Message] Fetching records", unit="records") as pbar:

        for batch in chunked(pmids, BATCH_SIZE):
            xml_raw  = fetch_xml(batch)
            records  = parse_pubmed_xml(xml_raw, next_pid)
            next_pid += len(records)

            for rec in records:
                json.dump(rec, out, ensure_ascii=False)
                out.write("\n")

            pbar.update(len(batch))

    print(f"[Message] Done. Wrote {next_pid-1} records to {OUT_FILE}")



if __name__ == "__main__":
    main()
