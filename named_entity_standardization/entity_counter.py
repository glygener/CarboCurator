import csv
import argparse
from collections import Counter

def count_frequencies(tsv_file):
    # Open the TSV file for reading.
    with open(tsv_file, 'r', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        # Exclude the "sentence_in_abstract" column from processing.
        #relevant_columns = [col for col in header if col != "sentence_in_abstract"]
        relevant_columns = [col for col in header if col == "Keywords"]
        # Create a counter for each relevant column.
        counters = {col: Counter() for col in relevant_columns}
        
        # Iterate over each row in the TSV file.
        for row in reader:
            # Skip rows that don't have the expected number of columns.
            if len(row) < len(header):
                continue
            # Process each relevant column.
            for col in relevant_columns:
                # Get the column index from the header.
                idx = header.index(col)
                # Ensure the index is valid for this row.
                if idx >= len(row):
                    continue
                value = row[idx]
                # Split the value by pipe and count each token.
                tokens = value.strip().split('|')
                for token in tokens:
                    token = token.strip()
                    if token:  # Only count non-empty tokens.
                        counters[col][token] += 1
    return counters

def main():
    parser = argparse.ArgumentParser(
        description='Count frequencies of values in each column of a TSV file (excluding "sentence_in_abstract").'
    )
    parser.add_argument('tsv_file', help='Path to the TSV file')
    args = parser.parse_args()

    frequencies = count_frequencies(args.tsv_file)

    # Write an output file for each column (excluding "sentence_in_abstract").
    for col, counter in frequencies.items():
        if col == "sentence_in_abstract":
            continue
        output_file = f"frequency_{col}.tsv"
        with open(output_file, 'w', encoding='utf-8') as outf:
            # Sort items by frequency in descending order.
            sorted_items = sorted(counter.items(), key=lambda x: x[1], reverse=True)
            outf.write(f"{col}\tfrequency\n")
            for token, count in sorted_items:
                outf.write(f"{token}\t{count}\n")

if __name__ == '__main__':
    main()