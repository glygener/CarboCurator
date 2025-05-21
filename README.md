# CarboCurator: Automated Glycan Biomarker Extraction

**Keywords**: `Text Mining`, `Glycobiology`, `Large Language Models`, `Python`

**Authors**: Cyrus C.H. Au Yeung  
**Date**: May 2025  
**Poster**: [poster.pdf](https://github.com/glygener/CarboCurator/blob/main/figures/poster.pdf#:~:text=poster.pdf)

> The project pipeline is being restructured to improve interoperability and incorporate new NER–RE and post-processing modules. Production-ready scripts now reside in `main/`, while beta versions of each major task live in their respective directories (e.g., `named_entity_recognition/`, `preprocessing/`). The standardization and mapping utilities for extracted glycan terms are housed in `named_entity_standardization/`.

---

## Project Overview

CarboCurator is a scalable pipeline for automated extraction and curation of **glycan–disease relationships** from biomedical literature. Leveraging PubMed abstracts and Large Language Models (LLMs), it recognizes glycan entities, associated diseases, specimen details, and evidence sentences, then structures the results according to a predefined JSON schema for downstream analysis and integration into biomarker knowledgebases.

**Objectives**

- Automate Named Entity Recognition (NER) of glycan entities, glycoproteins, glycoenzymes, diseases, specimen types, and biomarker roles.
- Link extracted entities through Relation Extraction (RE) to form structured glycan–disease associations.
- Map recognized terms to standard ontologies (GlyTouCan, Uberon, Disease Ontology) to enhance interoperability.
- Provide high‐quality, JSON‐formatted outputs for integration into knowledgebases.

---

## Data Source

- **PubMed Abstracts**: 4,595 articles retrieved via a custom Entrez query, including titles, abstracts, keywords, MeSH and chemical names.
- **Glycan Structure Dictionary** ([Vora *et al.*, 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10243773/)): for mapping glycan terms to standardized motifs.

---

## Methods

![workflow_diagram](https://github.com/glygener/CarboCurator/blob/main/figures/workflow_diagram.png#:~:text=workflow_diagram.png)

1. **Corpus Retrieval**  
   - Use Entrez E-utilities to fetch PubMed XML records.  
   - Filter out duplicates and entries lacking abstracts.

2. **Named Entity Recognition (NER)**  
   - Prompt‐driven extraction of five entity categories:  
     - **assessed_biomarker_entity** (glycan, glycan change, glycoprotein, glycoenzyme)  
     - **associated_disease** (disease name, medical intervention)  
     - **sample** (organism, specimen type)  
     - **biomarker_role** (diagnostic, prognostic, etc.)  
     - **evidence** (PMID + sentence)

3. **Relation Extraction (RE)**  
   - Link entities within sentences to form glycan–disease associations, capturing context and semantic roles.

4. **Standardization & Mapping**  
   - Map extracted names to GlyTouCan accessions, Uberon IDs, and Disease Ontology terms via MeSH/Chemical names.

5. **Output & Validation**  
   - Produce JSONL outputs conforming to the project schema.  
   - Convert to TSV for summary tables.

## Example

![example_curation](https://github.com/glygener/CarboCurator/blob/main/figures/curation_example.png#:~:text=curation_example.png)

---

## Project Structure

```bash
carbocurator/  
├── data/                        # All data files used in the project  
│   ├── raw/                     # Raw input data files  
│   │   ├── nen_dict.py          # Named Entity Normalization dictionary  
│   │   ├── ner_re/              # Raw NER and RE data  
│   │   └── nes/                 # Raw NES data  
│   └── processed/               # Processed data outputs  
│       ├── biomarker_table/     # Processed biomarker information  
│       ├── corpus/              # Processed text corpus  
│       ├── ner_re/              # NER & RE results  
│       └── nes/                 # Named Entity Standardization results  

├── main/                        # Core application scripts  
│   ├── abstracts_curate.py      # Abstract curation script  
│   ├── abstracts_curate.sh      # Shell script for curation pipeline  
│   ├── abstracts_validate.py    # Abstract validation script  
│   ├── json_parse.py            # JSON parsing utilities  
│   ├── pmid_fetch.py            # PubMed ID fetching script  
│   ├── preprocess_get_abst.py   # Abstract preprocessing  
│   ├── preprocess_get_pmid.py   # PMID preprocessing  
│   ├── requirements.txt         # Python dependencies for main  
│   ├── subset_pmids.txt         # Subset of PMIDs for testing  
│   └── temp_abstract.md         # Temporary abstract storage  

├── preprocessing/               # Data preprocessing scripts  
│   ├── preprocess_get_abst.py   # Abstract retrieval & processing  
│   ├── preprocess_get_pmid.py   # PMID retrieval & processing  
│   ├── preprocess_label_abst.py # Abstract labeling  
│   └── preprocess_qc_abst.py    # Abstract quality control  

├── named_entity_recognition/    # NER processing scripts  
│   ├── ner_bern2_qc.py          # BERN2 quality control  
│   ├── ner_bern2.py             # BERN2 processing  
│   ├── ner_biobert.py           # BioBERT processing  
│   ├── postprocess_bern2.py     # BERN2 post-processing  
│   ├── postprocess_biobert.py   # BioBERT post-processing  
│   └── tok_biobert.py           # BioBERT tokenization  

├── relation_extraction/         # Relation extraction processing  
│   ├── re_gpt.py                # GPT-based relation extraction  
│   ├── re_parse_json.py         # JSON parsing for relations  
│   └── re_summary.py            # Relation summary generation  

├── named_entity_standardization/ # Entity standardization scripts  
│   ├── entity_counter.py        # Entity frequency counter  
│   ├── gsd_extractor.py         # Gold Standard Data extractor  
│   └── term_matcher.py          # Term matching utilities  

├── model_validation/            # Validation scripts for models  
│   ├── ner_validation_auto.py   # Automated NER validation  
│   ├── ner_validation_manual.py # Manual NER validation  
│   ├── re_validation_auto.py    # Automated RE validation  
│   └── re_validation_manual.py  # Manual RE validation  

├── requirements.txt             # Project-wide Python dependencies  

└── README.md                    # Project overview and instructions  
```

## Major Dependencies

This project is written in **Python (≥ 3.8)** and uses:

- `openai` – for GPT‐based NER/RE  
- `biopython` – for Entrez API access  

Install with:

```bash
pip install -r requirements.txt
```
## How to Run

1. Clone the repo
```bash
git clone https://github.com/glygener/CarboCurator.git
cd CarboCurator
```
2. Set API keys ([NCBI](https://support.nlm.nih.gov/kbArticle/?pn=KA-05317)) ([OpenAI](https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://platform.openai.com/api-keys&ved=2ahUKEwicptvz97SNAxWsF1kFHecXN7cQFnoECBYQAQ&usg=AOvVaw1YhcGDWJXhiKSfmL59Pnfn))

Unix & macOS (bash/zsh):  
- Add to ~/.bashrc or ~/.zshrc
```bash
export NCBI_API_KEY="your_ncbi_api_key_here"
export OPENAI_API_KEY="your_openai_api_key_here"
```

- Then reload:
```bash
source ~/.bashrc
# or
source ~/.zshrc
```

Windows (PowerShell):  
```ps
[Environment]::SetEnvironmentVariable('NCBI_API_KEY','your_ncbi_api_key_here','User')
[Environment]::SetEnvironmentVariable('OPENAI_API_KEY','your_openai_api_key_here','User')
```
- Restart your terminal for changes to take effect

4. Prepare data
- Place your PubMed XML dump (or Entrez fetch script output) under `data/raw/`.
- Ensure `requirements.txt` is up to date and install dependencies.

5. Execute the pipeline
```bash
# Fetch and preprocess abstracts
python main/pmid_fetch.py
python main/preprocess_get_abst.py

# NER & RE
python main/abstracts_curate.py

# Postprocess
python main/json_parse.py
```

6. Inspect outputs
- JSONL files under `data/processed/ner_re/`
- Standardized TSV under `data/processed/nes/`

## License

This project is licensed under the [MIT License]().
