# the file number has to come right before '.json' in the input files, e.g. "model-4o-mini_output1.json"
import os
import json
import csv
import glob
import re
from optparse import OptionParser
from tqdm import tqdm

def main():
    #Option parser
    usage = "usage: %prog [-i input] [-o output]"
    parser = OptionParser(usage)
    parser.add_option("-i", action = "store", dest = "input",
        default = "extracted_abstracts",
        help = "Extracted abstracts directory name")
    parser.add_option("-o", action = "store", dest = "output",
        default = "biomarker_table.tsv",
        help = "Output biomarker table name")
    (options, args) = parser.parse_args()

    if options.output.endswith('.tsv') == False:
        options.output += '.tsv'

    # Paths
    script_dir = os.path.dirname(os.path.realpath(__file__))
    input_dir = os.path.normpath(os.path.join(script_dir, "../data_output" ,options.input))
    if not os.path.exists(input_dir):
        print(f"Error: '{os.path.basename(input_dir)}' not found.")
        parser.print_help()
        exit(0)
    json_files_path = input_dir + "/extracted_abstract*.json"
    os.makedirs(os.path.normpath(os.path.join(script_dir, "../data_output/biomarker_table")), exist_ok = True)
    output_path = os.path.normpath(os.path.join(script_dir, "../data_output/biomarker_table", options.output))
    if os.path.exists(output_path):
        print(f"Warning: '{os.path.basename(output_path)}' already exists. Do you want to overwrite it? (y/n)")
        response = input()
        if response.lower() != 'y' and response.lower() != 'yes':
            print("Exiting...")
            exit(0)
    
    json_list = glob.glob(json_files_path)
    sorted_json_list = sorted(json_list, key = lambda f: int(re.search(r'(\d+)(?=\.json$)', f).group())) # numberically sorted not lexicographically
    #print(sorted_json_list)

    tsv_columns = [
        "glycan_entity",
        "glycan_entity_change",
        "glycoprotein_entity_name",
        "glycoprotein_gene_name",
        "glycoprotein_mapped",
        "glycoenzyme_entity_name",
        "glycoenzyme_gene_name",
        "glycoenzyme_mapped",
        "biomarker_role",
        "disease_name",
        "disease_mapped",
        "medical_intervention",
        "organism",
        "specimen_type",
        "specimen_mapped",
        "pmid",
        "sentence_in_abstract"
    ]
    writetsv(sorted_json_list, output_path, tsv_columns)
    print(f"Biomarker table created: {os.path.basename(output_path)} in 'data_output/biomarker_table'")


def writetsv(sorted_json_list, output_path, tsv_columns):

    with open(output_path, 'w', newline = '', encoding = 'utf-8') as tsvfile:
        writer = csv.DictWriter(tsvfile, fieldnames = tsv_columns, delimiter = '\t')
        writer.writeheader()

        for json_file in tqdm(sorted_json_list, bar_format = "{l_bar}  |{bar} {n_fmt}/{total_fmt} [{elapsed}<{remaining}] {postfix}", desc = "Processing JSON files"):
            with open(json_file, 'r') as jf:
                data = json.load(jf)
                for biomarker in data.get('glycan_biomarkers', []):
                    row = {
                        "glycan_entity": biomarker["assessed_biomarker_entity"]["glycan_entity"],
                        "glycan_entity_change": biomarker["assessed_biomarker_entity"]["glycan_entity_change"],
                        "glycoprotein_entity_name": biomarker["assessed_biomarker_entity"]["glycoprotein_entity_name"],
                        "glycoprotein_gene_name": biomarker["assessed_biomarker_entity"]["glycoprotein_gene_name"],
                        "glycoprotein_mapped": biomarker["assessed_biomarker_entity"]["glycoprotein_mapped"],
                        "glycoenzyme_entity_name": biomarker["assessed_biomarker_entity"]["glycoenzyme_entity_name"],
                        "glycoenzyme_gene_name": biomarker["assessed_biomarker_entity"]["glycoenzyme_gene_name"],
                        "glycoenzyme_mapped": biomarker["assessed_biomarker_entity"]["glycoenzyme_mapped"],
                        "disease_name": biomarker["associated_disease"]["disease_name"],
                        "disease_mapped": biomarker["associated_disease"]["disease_mapped"],
                        "medical_intervention": biomarker["associated_disease"]["medical_intervention"],
                        "biomarker_role": "|".join(biomarker.get("biomarker_role", [])),
                        "organism": biomarker["sample"]["organism"],
                        "specimen_type": biomarker["sample"]["specimen_type"],
                        "specimen_mapped": biomarker["sample"]["specimen_mapped"],
                        "pmid": biomarker["evidence"]["pmid"],
                        "sentence_in_abstract": biomarker["evidence"]["sentence_in_abstract"]
                    }
                    writer.writerow(row)

if __name__ == "__main__":
    main()