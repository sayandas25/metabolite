# Metabolite Online Annotation & Literature Search App

This is the corrected Streamlit Cloud-ready version.

## Main correction

This app does not stop after checking the local annotation table. Even if the local annotation is missing, it still searches:

- PubChem
- KEGG
- ClassyFire
- PubMed
- Europe PMC

## Files needed for Streamlit Community Cloud

Upload these files to your GitHub repository root:

- app.py
- requirements.txt

Optional:

- metabolite_annotation_template.csv

## Deploy on Streamlit Community Cloud

1. Push the files to GitHub.
2. Go to https://share.streamlit.io/
3. Click New app.
4. Select your repository.
5. Set main file path as:

app.py

6. Deploy.

## Test example

Use:

Metabolite: Cysteine
Relation topic: preeclampsia

You should now see PubChem and KEGG fields populated, even if local annotation is not available.
