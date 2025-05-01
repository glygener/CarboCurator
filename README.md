# CarboCurator - Glycan-Disease Text Mining Pipeline

This repository contains a set of tools for mining glycan-related biomarker information from PubMed abstracts. The pipeline includes a PubMed abstract fetcher and a ChatGPT API-based NER tool to extract biological entities—specifically glycan epitopes/structural terms and disease associations.

## Features

- **PubMed Abstract Fetcher**  
  Fetches abstracts using a custom query and outputs a table with PMIDs, titles, abstracts, keywords, and MeSH IDs.

- **ChatGPT API Integration**  
  Processes abstracts with the ChatGPT API to identify and extract glycan-related entities and disease associations.

## Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/yourusername/your-repo.git
   cd your-repo

2. **Install dependencies:**
    ```bash
    
    ```

3. **Configure API keys:**
    ```bash
    export "CHATGPT_API_KEY"="{YOUR_API_KEY_HERE}"
    ```

## Usage
`query_builder.py`
`abstracts_fetch.py`
`abstracts_curate.py`
`json_parse.py`

## Project Structure

## Contributing

## License
