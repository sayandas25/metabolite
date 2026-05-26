# Metabolite Annotation v4 - No Blank Outputs

This is the corrected fresh version.

## Important
Use this version instead of all previous zip files.

## What it does
- Searches local curated table
- Searches PubChem with multiple name variants
- Searches KEGG with multiple name variants
- Searches ChEBI via EBI OLS
- Searches ClassyFire if InChIKey is available
- Searches PubMed and Europe PMC
- Avoids blank key outputs
- Gives closest likely organ/physiological compartment using fallback rules
- Classifies unresolved compounds as likely exogenous/xenobiotic or unclassified

## Streamlit Community Cloud deployment
Upload only these two required files to GitHub root:
- app.py
- requirements.txt

Then redeploy/reboot your Streamlit app.

## Test
Metabolite: Cysteine
Relation topic: preeclampsia

Expected:
Cysteine should now map to sulfur amino acid/redox metabolism and show a physiological compartment such as:
Blood/plasma; liver; cytosol; placenta; immune cells
