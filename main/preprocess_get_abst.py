import os
from Bio import Entrez
import pandas as pd
from optparse import OptionParser

def main():
    # Option parser
    usage = "usage: %prog [-i input] [-o output] [-b batchsize] [-e email]"
    parser = OptionParser(usage)
    parser.add_option("-i", action="store", dest="input",
        help="PMID list: input file name")
    parser.add_option("-o", action="store", dest="output",
        default="pubmed_abstracts",
        help="Abstracts: output directory name")
    parser.add_option("-b", action="store", dest="batch_size",
        type="int", default=1000,
        help="Batch size of PMIDs to be fetched [default: %default]")
    parser.add_option("-e", action="store", dest="email",
        default="your.username@example.com",
        help="Optional. Email address provided to NCBI E-utilities")
    (options, args) = parser.parse_args()

    if (not options.input):
        print("\033[91mError\033[0m: Input file must be provided.")
        parser.print_help()
        exit(0)
    
    # Paths
    script_dir = os.path.dirname(os.path.realpath(__file__))
    input_path = os.path.join(script_dir, "..", "data_input" ,options.input)
    input_path = os.path.normpath(input_path)
    output_dir = os.path.join(script_dir, "..", "data_output", options.output) # For storing abstract text files.
    output_dir = os.path.normpath(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(script_dir, "..", "data_output", options.output + ".tsv")
    output_path = os.path.normpath(output_path)

    Entrez.email = options.email # NCBI will attempt to contact a user at the email address provided before blocking access to the E-utilities.
    batch_size = options.batch_size # Number of PubMed articles fetched per batch.
    
    # Main
    try: 
        with open(input_path, 'r') as file:
            pmid_list = [pmid.strip() for pmid in file]
        print(f"[\033[92mProcess Started\033[0m] Batch size = {batch_size}. Fetching {len(pmid_list)} PMIDs... ")

    except FileNotFoundError:
        print(f"[\033[91mError\033[0m] File not found: {input_path}")
        parser.print_help()
        exit(0)

    tsv = []
    row = 1 # Index; row number in tsv data.

    for i in range(0, len(pmid_list), batch_size):
        batch_pmids = pmid_list[i:i+batch_size]
        xml = fetch_xml(batch_pmids)

        for article in xml["PubmedArticle"]:
            extracted_fields = extract_fields(article)
            save_txt(output_dir, extracted_fields, row)
            tsv.append(extracted_fields)
            row += 1
        
        batch_num = int(1+(i/batch_size))
        progress = min((i+batch_size)/len(pmid_list)*100, 100)
    
        print(f" Finish parsing batch {batch_num}.\tProgress: {progress:.1f}%")
    
    save_tsv(tsv, output_path)
    print(f"[\033[92mComplete\033[0m] Abstracts fetched and parsed to TSV and TXT.\n Text files saved to:\t{output_dir}\n TSV table saved to:\t{output_path}")

# Fetch data through Entrez in XML format
def fetch_xml(pmids):
    handle = Entrez.efetch(db="pubmed", id=",".join(pmids), retmode="xml")
    xml = Entrez.read(handle)
    handle.close()
    return xml

# Extract fields of interest from XML data
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
    df.to_csv(output_path, sep = '\t', index = False)


def save_txt(output_dir, extracted_fields, row):
    abs_filename = os.path.join(output_dir, f"pubmed_abstract{row}.txt") # abstracts_curate.py dependency
    with open(abs_filename, "w", encoding = "utf-8") as abs_file: # abs means abstract
        abs_file.write(
            f"PMID: {extracted_fields['PMID']}\n"
            f"Title: {extracted_fields['Title']}\n"
            f"Abstract: {extracted_fields['Abstract']}\n"
            f"Keywords: {extracted_fields['Keywords']}\n"
            f"Chemical names: {extracted_fields['Chemical_names']}\n"
            f"MeSH names: {extracted_fields['MeSH_names']}\n"
        )


if __name__ == "__main__":
    main()
