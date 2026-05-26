# Metabolite Annotation & Literature Search Streamlit App

## Purpose
This app allows the user to enter one metabolite or a list of metabolites and retrieve:

1. Compound class
2. Pathway
3. Superpathway
4. Physiological compartment
5. PubMed papers linking the metabolite to a user-defined disease, exposure, phenotype, or outcome

## How to run locally

1. Install Python 3.10 or later.
2. Open Command Prompt or Terminal.
3. Go to this folder.
4. Run:

pip install -r requirements.txt

5. Run:

streamlit run app.py

## Input options
You can use the small demo annotation database included in the app, or upload your own CSV.

Your CSV should contain these columns:

metabolite, synonyms, compound_class, pathway, superpathway, physiological_compartment, hmdb_id, kegg_id, notes

## Recommended future improvements
- Add HMDB XML parsing
- Add KEGG REST API
- Add ClassyFire annotation by InChIKey/SMILES
- Add LitSense/Semantic Scholar abstracts
- Add evidence grading: human study, animal study, in vitro, review
- Add AI-generated biological interpretation with citations
