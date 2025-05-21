# Script 3: abstracts_curate.py
# Description: This script curates abstracts from PubMed using OpenAI's API. It extracts relevant information based on a provided JSON schema and query text.
# Usage: python abstracts_curate.py -i <input_dir> -o <output_dir> -s <schema_file> -q <query_file> -k <api_key>
# Input: Abstract text files in the data_input directory, a JSON schema, and a query text file.
# Output: JSON files containing the extracted entities in the data_output directory.
# Dependencies: openai, tqdm
# Note: Make sure to set the OpenAI API key as an environment variable or provide it through command line arguments.
#       - To set the OpenAI API key as an environment variable, use:
#         $ export 'OPENAI_API_KEY'='sk...'
#       - To automatically set the OpenAI API key as environmental variable, store it in ~/.bashrc (linux) or ~/.zshrc (macOS)
#       - To provide the OpenAI API key through command line arguments, use:
#         $ python abstracts_curate.py -k 'sk...'
#       Make sure the data_input and data_output directories exist before running, and make sure both the json schema and query text file are both at data_input.
#       User may use run_curate.sh to run process in the background. Realtime logs file available.
#       For rate limit errors, reduce the number of workers in parallel processing (max_workers).

__version__ = "1.3"
__release__ = "2025-03-31"

import os
import json
import time
import glob
import re
import optparse
import logging
import concurrent.futures
from itertools import repeat

import openai
from openai import OpenAI
from tqdm import tqdm

###############################################################################
# Existing helper -------------------------------------------------------------
###############################################################################

def import_api_key(options):
    if not options.api_key:
        api_key = os.getenv("OPENAI_API_KEY")
        if not api_key:
            print("API key not found. Please provide API key using -k flag or set as environmental variable.")
            parser.print_help()
            exit(0)
    else:
        api_key = options.api_key
    return api_key

# GPT models. To see all available models, visit https://platform.openai.com/docs/models
model_list = {
    "4.1": "gpt-4.1-2025-04-14", "4.1-mini": "gpt-4.1-mini-2025-04-14", "4.1-nano": "gpt-4.1-nano-2025-04-14",
    "4o": "gpt-4o-2024-08-06", "4o-mini": "gpt-4o-mini-2024-07-18",
    "o3": "o3-2025-04-16", "o3-mini": "o3-mini-2025-01-31",
    "o4-mini": "o4-mini-2025-04-16"
}  # List updated on 2025-05-02

# -----------------------------------------------------------------------------
# The following OpenAI helper is intentionally left untouched beyond schema fix
# -----------------------------------------------------------------------------

def openai_api(api_key):
    response = OpenAI(api_key=api_key).responses.create(
        model=model_list["4.1"],
        store=False,
        temperature=0,
        input=[
            {"role": "system", "content": "You are a biocurator expert in named entity recognition for glycobiology and biomedical literature."},
            {"role": "user",   "content": prompt}
        ],
        text={
            "format": {
                "type": "json_schema",
                "name": "text_mining",
                "strict": True,
                "schema": schema_obj
            }
        }
    )
    print(response.usage.total_tokens)
    return response.output_text

###############################################################################
# Modified ner_re -------------------------------------------------------------
###############################################################################

def ner_re(api_key, query, abstract_text, schema_obj):
    globals()["prompt"] = f"{query}\nAbstract text:\n{abstract_text}"
    globals()["schema_obj"] = schema_obj

    max_retries = 5
    retries_rate = 0
    retries_conn = 0

    while True:
        try:
            response_text = openai_api(api_key)
            try:
                return json.loads(response_text)
            except json.JSONDecodeError as e:
                with open("error.json", "w") as f:
                    f.write(response_text)
                print(f"JSON Decode Error: {e}")
                raise Exception("JSON decode error")
        except openai.RateLimitError as e:
            if retries_rate < max_retries:
                retries_rate += 1
                print(f"Rate limit exceeded. Retrying in 30 seconds... ({retries_rate}/{max_retries})")
                time.sleep(30)
            else:
                print("Rate limit exceeded. Maximum retries reached. Exiting.")
                exit(0)
        except openai.APIConnectionError as e:
            if retries_conn < max_retries:
                retries_conn += 1
                print(f"Connection error. Retrying in 30 seconds... ({retries_conn}/{max_retries})")
                time.sleep(30)
            else:
                print("Connection error. Maximum retries reached. Exiting.")
                exit(0)

###############################################################################
# CLI / path helpers ----------------------------------------------------------
###############################################################################

def parse_options():
    global parser
    usage = "usage: %prog [-i input] [-o output] [-s schema] [-q query] [-k api_key] [-w workers]"
    parser = optparse.OptionParser(usage)

    parser.add_option("-i", action="store", dest="input", default="abstracts", help="abstract text directory name")
    parser.add_option("-o", action="store", dest="output", default="extracted_abstracts", help="Output directory path")
    parser.add_option("-s", action="store", dest="schema", default="schema-v5.0_test.json", help="JSON schema path")
    parser.add_option("-q", action="store", dest="query", default="query-v3.0_test.txt", help="Query text path")

    group = optparse.OptionGroup(parser, "API Key Options",
                                 "You may supply the OPENAI API key with -k or as env var OPENAI_API_KEY")
    group.add_option("-k", action="store", dest="api_key", help="OpenAI API key")
    parser.add_option_group(group)

    parser.add_option("-w", action="store", dest="workers", type="int", default=15, help="Number of parallel workers (default 15)")
    return parser.parse_args()[0]


def build_paths(options):
    script_dir = os.path.dirname(os.path.realpath(__file__))
    input_dir = os.path.normpath(os.path.join(script_dir, "..", "data_output", options.input))
    schema_path = os.path.normpath(os.path.join(script_dir, "..", "data_input", options.schema))
    query_path = os.path.normpath(os.path.join(script_dir, "..", "data_input", options.query))
    output_dir = os.path.normpath(os.path.join(script_dir, "..", "data_output", options.output))
    return input_dir, output_dir, schema_path, query_path


def list_index_strings(input_dir):
    files = glob.glob(os.path.join(input_dir, "abstract_??????.txt"))
    idx_strings = []
    for f in files:
        m = re.search(r"abstract_(\d{6})\.txt$", os.path.basename(f))
        if m:
            idx_strings.append(m.group(1))
    if not idx_strings:
        raise FileNotFoundError("No abstract files found.")
    return sorted(idx_strings)

###############################################################################
# Abstract processor ---------------------------------------------------------
###############################################################################

def process_single_abstract(idx_str, paths, api_key, query, schema_obj, log):
    start_ts = time.time()
    input_dir, output_dir = paths
    input_file = os.path.join(input_dir, f"abstract_{idx_str}.txt")
    output_file = os.path.join(output_dir, f"extracted_{idx_str}.json")

    if not os.path.exists(input_file):
        log.warning(f"{time.strftime('%Y-%m-%d %H:%M:%S')}\t{idx_str}\tMISSING_INPUT")
        return

    with open(input_file, "r") as fh:
        abstract_text = fh.read()

    m_ab = re.search(r"^Abstract:\s*(.*)", abstract_text, flags=re.MULTILINE | re.DOTALL)
    if m_ab and not m_ab.group(1).strip():
        with open(output_file, "w") as jf:
            json.dump({}, jf)
        elapsed = time.time() - start_ts
        log.info(f"{time.strftime('%Y-%m-%d %H:%M:%S')}\t{idx_str}\tSKIPPED_EMPTY\t{elapsed:.1f}s")
        return

    result = ner_re(api_key, query, abstract_text, schema_obj)

    with open(output_file, "w") as jf:
        json.dump(result, jf, indent=2)

    elapsed = time.time() - start_ts
    log.info(f"{time.strftime('%Y-%m-%d %H:%M:%S')}\t{idx_str}\tPROCESSED\t{elapsed:.1f}s")

###############################################################################
# Main -----------------------------------------------------------------------
###############################################################################

def main():
    options = parse_options()
    api_key = import_api_key(options)
    input_dir, output_dir, schema_path, query_path = build_paths(options)

    for name, pth in [("input directory", input_dir), ("schema", schema_path), ("query", query_path)]:
        if not os.path.exists(pth):
            print(f"{name.capitalize()} '{pth}' does not exist.")
            exit(0)

    os.makedirs(output_dir, exist_ok=True)

    with open(schema_path, "r") as fh:
        schema_obj = json.load(fh)
    with open(query_path, "r") as fh:
        query = fh.read()

    indices = list_index_strings(input_dir)

    log_path = os.path.join(output_dir, "curation.log")
    logging.basicConfig(filename=log_path, filemode="a", format="%(message)s", level=logging.INFO)
    log = logging.getLogger(__name__)

    log.info("#-----------------------------------------------------------")
    log.info(f"[Process Started]  Time: {time.ctime()}")
    log.info(f" Abstract directory: {os.path.basename(input_dir)}")
    log.info(f" Output directory:   {os.path.basename(output_dir)}")
    log.info(f" Schema:             {os.path.basename(schema_path)}")
    log.info(f" Query:              {os.path.basename(query_path)}")
    log.info("#-----------------------------------------------------------")

    start_all = time.time()
    paths = (input_dir, output_dir)

    with concurrent.futures.ThreadPoolExecutor(max_workers=options.workers) as executor:
        iterator = executor.map(process_single_abstract,
                                indices,
                                repeat(paths),
                                repeat(api_key),
                                repeat(query),
                                repeat(schema_obj),
                                repeat(log))
        for _ in tqdm(iterator, total=len(indices), desc="Curating text", colour="green", unit="abstract"):
            pass

    elapsed_all = (time.time() - start_all) / 60.0
    log.info(f"Finished batch at {time.strftime('%Y-%m-%d %H:%M:%S')} – total elapsed {elapsed_all:.1f} min")


if __name__ == "__main__":
    main()
