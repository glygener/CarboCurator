# Make sure the data_input and data_output directories exist before running and make sure both the json schema and query text file are both at data_input
# Use run_curate.sh to run process in background. Realtime log file available.
import os
import json
import time
import glob, re
import optparse
import openai
import concurrent.futures
from openai import OpenAI
from tqdm import tqdm

__version__ = "1.3"
__release__ = "2025-03-31"

def main():
    # Option parser
    usage = "usage: %prog [-i input] [-o output] [-s schema] [-q query] [-k api_key]"
    parser = optparse.OptionParser(usage)
    parser.add_option("-i", action="store", dest="input",
        default="pubmed_abstracts",
        help="abstract text directory name")
    parser.add_option("-o", action="store", dest="output",
        default="extracted_abstracts",
        help="Output directory path")
    parser.add_option("-s", action="store", dest="schema",
        default="schema-v5.0_test.json",
        help="JSON schema path")
    parser.add_option("-q", action="store", dest="query",
        default="query-v3.0_test.txt",
        help="Query text path")

    group = optparse.OptionGroup(parser, "API Key Options",
        "You may choose to supply the ChatGPT API key through the -k flag, "
        "or load key as environmental variable: "
        "$ export 'CHATGPT_API_KEY'='sk...'")
    group.add_option("-k", action="store", dest="api_key",
        help="ChatGPT API key")
    parser.add_option_group(group)

    (options, args) = parser.parse_args()

    # Paths
    script_dir = os.path.dirname(os.path.realpath(__file__))
    input_dir = os.path.normpath(os.path.join(script_dir, "..", "data_output", options.input))
    schema_path = os.path.normpath(os.path.join(script_dir, "..", "data_input", options.schema))
    query_path = os.path.normpath(os.path.join(script_dir, "..", "data_input", options.query))
    output_dir = os.path.normpath(os.path.join(script_dir, "..", "data_output", options.output))

    # Check if API key is provided
    if not options.api_key:
        api_key = os.getenv("CHATGPT_API_KEY")
        if not api_key:
            print("API key not found. Please provide API key using -k flag or set as environmental variable.")
            parser.print_help()
            exit(0)
    else:
        api_key = options.api_key

    # Check if input, schema, and query files exist
    for opt, opt_path in [("input", input_dir), ("schema", schema_path), ("query", query_path)]:
        if not os.path.exists(opt_path):
            print(f"{opt.capitalize()} file '{opt_path}' does not exist.")
            parser.print_help()
            exit(0)

    os.makedirs(output_dir, exist_ok=True) # Do not reposition this line

    # Find largest abstract index
    input_files_path = os.path.join(input_dir, "pubmed_abstract*.txt")
    files = glob.glob(input_files_path)
    indices = []
    for file in files:
        match = re.search(r'(\d+)\.txt', os.path.basename(file))
        if match:
            indices.append(int(match.group(1)))
    if not indices:
        print("No abstract files found. File names should be formatted as 'pubmed_abstract*.txt'")
        exit(0)
    max_idx = max(indices)

    # Load schema and query
    with open(schema_path, "r") as schema_file:
        schema_obj = json.load(schema_file)
    with open(query_path, "r") as query_file:
        query = query_file.read()

    start = time.time()
    print("[Process Started]")
    print(f" Time:               {time.ctime()}")
    print(f" Abstract directory: {os.path.basename(input_dir)}")
    print(f" Output directory:   {os.path.basename(output_dir)}")
    print(f" Schema:             {os.path.basename(schema_path)}")
    print(f" Query:              {os.path.basename(query_path)}")
    print("-----------------------------------------------------------")
    ##############################################################################
    def parallel_extraction(i): # nested within api_call()
        input_abs_path = os.path.join(input_dir, f"pubmed_abstract{i}.txt")
        output_json_path = os.path.join(output_dir, f"extracted_abstract{i}.json")

        if not os.path.exists(input_abs_path):
            return f"Omitting abstract {i}. File not found: {input_abs_path}"

        with open(input_abs_path, "r") as file:
            abstract_text = file.read()

        result_json = api_call(api_key, query, abstract_text, schema_obj)

        with open(output_json_path, "w") as json_file:
            json.dump(result_json, json_file, indent=2)
    ##############################################################################

    start_idx = 1
    max_idx = 100 #########################REMOVE LINE 
    total_tasks = max_idx - start_idx + 1
    max_workers = 15 ####################ADJUST IF NEEDED (DEFAULT = 15). RATE IS LIMITED BY USAGE TIER. See https://platform.openai.com/settings/organization/limits

    with concurrent.futures.ThreadPoolExecutor(max_workers = max_workers) as executor:
        iterator = executor.map(parallel_extraction, range(start_idx, max_idx + 1))
        for result in tqdm(iterator, total = total_tasks, desc = "Processing abstracts", colour = "green"):
            pass
    end = time.time()
    print(f"Abstract curation completed. Elapsed time: {((end - start) / 60):.01f} min")

def api_call(api_key, query, abstract_text, schema_obj):
    prompt = f"{query}\nAbstract text:\n{abstract_text}"
    client = OpenAI(api_key=api_key)
    max_retries = 5  
    retries_rate = 0 # rate limit
    retries_conn = 0 # connection error
    # 'temperature' is not supported with the o3-mini model
    model_list = {"4o":"gpt-4o-2024-11-20", "4o-FT":"gpt-4o-2024-08-06","4o-mini":"gpt-4o-mini-2024-07-18", "o3-mini":"o3-mini-2025-01-31"}
    while True:
        try:
            response = client.chat.completions.create(
                model = model_list["o3-mini"],  # Use the appropriate model version
                store = False,
                #metadata = {"description": "biomarker_20250331_4o"},
                #temperature = 0,
                reasoning_effort = "medium",
                messages = [
                    {"role": "system", "content": "You are a biocurator expert in named entity recognition for glycobiology and biomedical literature."},
                    {"role": "user", "content": prompt}
                ],
                response_format = {
                    "type": "json_schema",
                    "json_schema": {
                        "name": "text_mining",
                        "strict": True,
                        "schema": schema_obj
                    }
                }
            )
            try:
                result_json = json.loads(response.choices[0].message.content)
                return result_json
            except json.JSONDecodeError as e:
                with open("error.json", "w") as f:
                    f.write(response.choices[0].message.content)
                print(f"JSON Decode Error: {e}")
                print(response.choices[0].message.content)
                print("Output length exceeds limit: check you query.\n")
                raise Exception("JSON decode error")
                exit(0) # Fix required
        except openai.RateLimitError as e:
            if retries_rate < max_retries:
                retries_rate += 1
                print(f"Rate limit exceeded. Retrying in 30 seconds... Retries: {retries_rate}/{max_retries}")
                time.sleep(30)
            else:
                print("Rate limit exceeded. Maximum retries reached. Exiting.")
                print(f"Error: {e}\n If problem persists, reduce the number of workers in parallel processing (max_workers).")
                exit(0)
        except openai.APIConnectionError as e:
            if retries_conn < max_retries:
                retries_conn += 1
                print(f"Connection error. Retrying in 30 seconds... Retries: {retries_conn}/{max_retries}")
                time.sleep(30)
            else:
                print("Connection error. Maximum retries reached. Exiting.")
                print(f"Error: {e}")
                exit(0)

if __name__ == "__main__":
    main()