from Bio import Entrez
import time

Entrez.email = "your.email@example.com"  # Required by NCBI

query = '((glycan*[tiab] OR oligosaccharide*[tiab] OR polysaccharide*[tiab] OR "glycoepitope"[tiab] OR "glycotopes"[tiab] OR "glycan motif"[tiab] OR "glycan substructure"[tiab] OR "glycan epitope"[tiab] OR "glycan determinant"[tiab] OR "glycan signature"[tiab] OR "glycan sequence"[tiab] OR "sialyl Lewis"[tiab] OR "Tn antigen"[tiab] OR "Thomsen-Friedenreich"[tiab] OR "Lewis X"[tiab] OR "Lewis A"[tiab] OR "H type 1"[tiab] OR "H type 2"[tiab] OR "blood group antigen"[tiab] OR "glycan code"[tiab] OR "glycan map"[tiab] OR "glycan fingerprint"[tiab])) AND hasabstract[text] AND ("journal article"[Publication Type] NOT review[Publication Type])'

retmax = 1
all_pmids = []

# Initial request to get total number of results
handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
record = Entrez.read(handle)
total = int(record["Count"])

# Loop over chunks of retmax
for start in range(0, total, retmax):
    handle = Entrez.esearch(db="pubmed", term=query, retstart=start, retmax=retmax)
    record = Entrez.read(handle)
    all_pmids.extend(record["IdList"])
    print(f"Fetched PMIDs {start + 1} to {start + len(record['IdList'])}")
    time.sleep(0.5)  # Respect rate limits

# Save PMIDs to file
with open("pmid_expanded.txt", "w") as f:
    for pmid in all_pmids:
        f.write(pmid + "\n")
