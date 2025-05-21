import csv
import json


def main():
    # Input and output file paths
    input_file = './data_input/gsd_v2.7.1_cleaned.tsv'
    output_file = './gsd_v2.7.1_cleaned.json'

    gsd_list = []
    # Read the TSV file
    with open(input_file, encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            ac = row.get('glytoucan_ac', '').strip()
            term = row.get('term', '').strip()
            syns = row.get('synonyms', '').strip()
            # Build synonyms list, always including the representative term
            if syns:
                synonyms = [term] + syns.split('|')
            else:
                synonyms = [term]
            # Create entry
            entry = {
                'glytoucan_ac': ac,
                'term': term,
                'synonyms': synonyms
            }
            gsd_list.append(entry)

    # Write out to JSON without nesting under a key
    with open(output_file, 'w', encoding='utf-8') as out:
        json.dump(gsd_list, out, indent=2, ensure_ascii=False)

if __name__ == '__main__':
    main()

