# Metabolite Annotation v5 - Gene and Protein Associations

This version adds a separate Gene/protein associations tab.

## New features
- Metabolite annotation with no blank key outputs
- Gene and protein association table
- KEGG EC enzyme number to UniProt human protein lookup
- UniProt enrichment for gene details
- Local curated gene/protein fallback
- Rule-based gene/protein fallback when no database association is found
- PubMed and Europe PMC literature search

## Required files for Streamlit Community Cloud
Upload these two files to your GitHub repository root:
- app.py
- requirements.txt

Optional:
- README_DEPLOY_V5.txt
- metabolite_annotation_template_v5.csv

## Suggested test
Metabolite:
Cysteine

Relation topic:
preeclampsia

Expected:
The annotation tab should show cysteine as a sulfur-containing amino acid/redox metabolite.
The gene/protein tab should include genes such as CBS, CTH, GCLC, GCLM, GSS, SLC7A11 and GGT1.
