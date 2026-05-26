
import streamlit as st
import pandas as pd
import requests
import xml.etree.ElementTree as ET
from io import StringIO
from urllib.parse import quote_plus
import time

# ----------------------------
# Page configuration
# ----------------------------
st.set_page_config(
    page_title="Metabolite Annotation & Literature Search App",
    page_icon="🧬",
    layout="wide"
)

# ----------------------------
# Helper functions
# ----------------------------
@st.cache_data(show_spinner=False)
def load_annotation_database(uploaded_file=None):
    """
    Loads metabolite annotation database.
    Expected columns:
    metabolite, synonyms, compound_class, pathway, superpathway,
    physiological_compartment, hmdb_id, kegg_id, notes
    """
    if uploaded_file is not None:
        return pd.read_csv(uploaded_file)

    # Small demo database. Replace or expand using HMDB/KEGG/LipidMaps exports.
    demo = pd.DataFrame([
        {
            "metabolite": "L-Arginine",
            "synonyms": "Arginine; L Arg; Arg",
            "compound_class": "Amino acid",
            "pathway": "Arginine and proline metabolism; Nitric oxide biosynthesis",
            "superpathway": "Amino acid metabolism",
            "physiological_compartment": "Plasma; cytosol; placenta; endothelial cells",
            "hmdb_id": "HMDB0000517",
            "kegg_id": "C00062",
            "notes": "Substrate for nitric oxide synthase and arginase."
        },
        {
            "metabolite": "L-Carnitine",
            "synonyms": "Carnitine; Levocarnitine",
            "compound_class": "Quaternary ammonium compound",
            "pathway": "Carnitine shuttle; fatty acid beta-oxidation",
            "superpathway": "Lipid and energy metabolism",
            "physiological_compartment": "Plasma; mitochondria; muscle; liver; placenta",
            "hmdb_id": "HMDB0000062",
            "kegg_id": "C00318",
            "notes": "Transports long-chain fatty acids into mitochondria."
        },
        {
            "metabolite": "Tryptophan",
            "synonyms": "L-Tryptophan; Trp",
            "compound_class": "Aromatic amino acid",
            "pathway": "Kynurenine pathway; serotonin biosynthesis; indole metabolism",
            "superpathway": "Amino acid metabolism; immune-metabolic regulation",
            "physiological_compartment": "Plasma; gut; placenta; immune cells; brain",
            "hmdb_id": "HMDB0000929",
            "kegg_id": "C00078",
            "notes": "Precursor of kynurenine, serotonin, melatonin, and microbial indoles."
        },
        {
            "metabolite": "Choline",
            "synonyms": "2-Hydroxyethyltrimethylammonium; Bilineurine",
            "compound_class": "Quaternary ammonium compound",
            "pathway": "Phosphatidylcholine metabolism; one-carbon metabolism",
            "superpathway": "Lipid metabolism; methyl donor metabolism",
            "physiological_compartment": "Plasma; liver; placenta; cell membrane",
            "hmdb_id": "HMDB0000097",
            "kegg_id": "C00114",
            "notes": "Important for membrane phospholipids and methylation biology."
        },
        {
            "metabolite": "Sphingosine-1-phosphate",
            "synonyms": "S1P; Sphingosine 1 phosphate",
            "compound_class": "Sphingolipid",
            "pathway": "Sphingolipid metabolism; vascular signaling",
            "superpathway": "Lipid signaling",
            "physiological_compartment": "Plasma; endothelial cells; immune cells; placenta",
            "hmdb_id": "HMDB0000277",
            "kegg_id": "C06124",
            "notes": "Bioactive lipid involved in vascular tone, angiogenesis, and immune trafficking."
        }
    ])
    return demo


def normalize_text(x):
    if pd.isna(x):
        return ""
    return str(x).strip().lower()


def find_metabolite_annotation(metabolite_name, db):
    """
    Searches metabolite name and synonyms.
    Returns best matching row or None.
    """
    query = normalize_text(metabolite_name)
    if query == "":
        return None

    db_copy = db.copy()
    db_copy["metabolite_norm"] = db_copy["metabolite"].apply(normalize_text)
    db_copy["synonyms_norm"] = db_copy.get("synonyms", "").apply(normalize_text)

    exact = db_copy[db_copy["metabolite_norm"] == query]
    if not exact.empty:
        return exact.iloc[0].to_dict()

    synonym_match = db_copy[db_copy["synonyms_norm"].str.contains(query, regex=False, na=False)]
    if not synonym_match.empty:
        return synonym_match.iloc[0].to_dict()

    partial = db_copy[
        db_copy["metabolite_norm"].str.contains(query, regex=False, na=False) |
        db_copy["synonyms_norm"].str.contains(query, regex=False, na=False)
    ]
    if not partial.empty:
        return partial.iloc[0].to_dict()

    return None


def pubmed_search(metabolite, relation_topic, max_results=10, email=None, api_key=None):
    """
    Searches PubMed using NCBI E-utilities.
    Returns a list of publication dictionaries.
    """
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    query = f'("{metabolite}"[Title/Abstract]) AND ("{relation_topic}"[Title/Abstract])'

    params = {
        "db": "pubmed",
        "term": query,
        "retmode": "json",
        "retmax": max_results,
        "sort": "relevance"
    }

    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key

    try:
        search_response = requests.get(base_url + "esearch.fcgi", params=params, timeout=20)
        search_response.raise_for_status()
        search_json = search_response.json()
        pmids = search_json.get("esearchresult", {}).get("idlist", [])
    except Exception as e:
        return [{"error": f"PubMed search failed: {e}"}]

    if not pmids:
        return []

    fetch_params = {
        "db": "pubmed",
        "id": ",".join(pmids),
        "retmode": "xml"
    }
    if email:
        fetch_params["email"] = email
    if api_key:
        fetch_params["api_key"] = api_key

    try:
        fetch_response = requests.get(base_url + "efetch.fcgi", params=fetch_params, timeout=30)
        fetch_response.raise_for_status()
        root = ET.fromstring(fetch_response.content)
    except Exception as e:
        return [{"error": f"PubMed fetch failed: {e}"}]

    papers = []
    for article in root.findall(".//PubmedArticle"):
        pmid = article.findtext(".//PMID", default="")
        title = article.findtext(".//ArticleTitle", default="")
        journal = article.findtext(".//Journal/Title", default="")
        year = article.findtext(".//PubDate/Year", default="")
        if not year:
            year = article.findtext(".//PubDate/MedlineDate", default="")

        abstract_parts = []
        for abstract_text in article.findall(".//Abstract/AbstractText"):
            label = abstract_text.attrib.get("Label", "")
            text = "".join(abstract_text.itertext())
            if label:
                abstract_parts.append(f"{label}: {text}")
            else:
                abstract_parts.append(text)

        abstract = " ".join(abstract_parts)

        authors = []
        for author in article.findall(".//Author")[:5]:
            last = author.findtext("LastName", default="")
            fore = author.findtext("ForeName", default="")
            name = f"{fore} {last}".strip()
            if name:
                authors.append(name)

        papers.append({
            "PMID": pmid,
            "Title": title,
            "Authors": "; ".join(authors),
            "Journal": journal,
            "Year": year,
            "Abstract": abstract,
            "PubMed_URL": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else ""
        })

    return papers


def make_interpretation(row, relation_topic, papers):
    if row is None:
        return "No local annotation was found. Consider adding this metabolite to the annotation database with HMDB/KEGG/LipidMaps identifiers."

    pathway = row.get("pathway", "")
    superpathway = row.get("superpathway", "")
    compartment = row.get("physiological_compartment", "")

    paper_count = len([p for p in papers if "error" not in p])
    return (
        f"This metabolite belongs mainly to {row.get('compound_class', 'an unspecified class')}. "
        f"It maps to {pathway} under the broader {superpathway}. "
        f"The reported physiological compartment/source includes {compartment}. "
        f"For the user-defined relation topic '{relation_topic}', the PubMed search returned {paper_count} article(s). "
        f"Manual biological review is recommended before using it as a biomarker claim."
    )


# ----------------------------
# Sidebar
# ----------------------------
st.sidebar.title("Settings")

uploaded_db = st.sidebar.file_uploader(
    "Upload metabolite annotation CSV",
    type=["csv"],
    help="CSV should contain: metabolite, synonyms, compound_class, pathway, superpathway, physiological_compartment, hmdb_id, kegg_id"
)

max_papers = st.sidebar.slider("Maximum PubMed papers per metabolite", 3, 30, 10)

email = st.sidebar.text_input(
    "NCBI email, optional but recommended",
    placeholder="your.email@example.com"
)

api_key = st.sidebar.text_input(
    "NCBI API key, optional",
    type="password"
)

annotation_db = load_annotation_database(uploaded_db)

# ----------------------------
# Main app
# ----------------------------
st.title("🧬 Metabolite Annotation & Literature Search App")

st.markdown("""
Enter one metabolite or paste a list of metabolites.  
The app returns **compound class, pathway, superpathway, physiological compartment**, and searches PubMed for papers linking each metabolite with your selected biological condition, disease, exposure, or outcome.
""")

col1, col2 = st.columns([1, 1])

with col1:
    relation_topic = st.text_input(
        "Search relation to:",
        placeholder="Example: preeclampsia, gestational hypertension, preterm birth, fetal growth restriction"
    )

with col2:
    input_mode = st.radio(
        "Input mode",
        ["Single metabolite", "List of metabolites"],
        horizontal=True
    )

if input_mode == "Single metabolite":
    metabolite_text = st.text_input(
        "Metabolite name",
        placeholder="Example: L-Arginine"
    )
else:
    metabolite_text = st.text_area(
        "Paste metabolite list",
        placeholder="One metabolite per line, or comma-separated list"
    )

run_button = st.button("Run annotation and literature search", type="primary")

if run_button:
    if not relation_topic.strip():
        st.error("Please enter the disease/outcome/condition in the 'Search relation to' box.")
        st.stop()

    if not metabolite_text.strip():
        st.error("Please enter at least one metabolite.")
        st.stop()

    # Parse metabolites
    raw_items = metabolite_text.replace(",", "\n").split("\n")
    metabolites = [x.strip() for x in raw_items if x.strip()]
    metabolites = list(dict.fromkeys(metabolites))  # remove duplicates while preserving order

    all_annotation_rows = []
    all_publication_rows = []

    progress = st.progress(0)
    status = st.empty()

    for i, metabolite in enumerate(metabolites):
        status.write(f"Processing: {metabolite}")

        annotation = find_metabolite_annotation(metabolite, annotation_db)
        papers = pubmed_search(
            metabolite=metabolite,
            relation_topic=relation_topic,
            max_results=max_papers,
            email=email,
            api_key=api_key
        )

        interpretation = make_interpretation(annotation, relation_topic, papers)

        if annotation is None:
            annotation_row = {
                "input_metabolite": metabolite,
                "matched_metabolite": "",
                "compound_class": "",
                "pathway": "",
                "superpathway": "",
                "physiological_compartment": "",
                "hmdb_id": "",
                "kegg_id": "",
                "notes": "",
                "interpretation": interpretation
            }
        else:
            annotation_row = {
                "input_metabolite": metabolite,
                "matched_metabolite": annotation.get("metabolite", ""),
                "compound_class": annotation.get("compound_class", ""),
                "pathway": annotation.get("pathway", ""),
                "superpathway": annotation.get("superpathway", ""),
                "physiological_compartment": annotation.get("physiological_compartment", ""),
                "hmdb_id": annotation.get("hmdb_id", ""),
                "kegg_id": annotation.get("kegg_id", ""),
                "notes": annotation.get("notes", ""),
                "interpretation": interpretation
            }

        all_annotation_rows.append(annotation_row)

        for p in papers:
            if "error" in p:
                all_publication_rows.append({
                    "metabolite": metabolite,
                    "error": p["error"]
                })
            else:
                p["metabolite"] = metabolite
                all_publication_rows.append(p)

        time.sleep(0.35)
        progress.progress((i + 1) / len(metabolites))

    status.write("Completed.")

    annotation_results = pd.DataFrame(all_annotation_rows)
    publication_results = pd.DataFrame(all_publication_rows)

    st.subheader("1. Metabolite annotation results")
    st.dataframe(annotation_results, use_container_width=True)

    st.download_button(
        "Download annotation results CSV",
        data=annotation_results.to_csv(index=False).encode("utf-8"),
        file_name="metabolite_annotation_results.csv",
        mime="text/csv"
    )

    st.subheader("2. PubMed literature search results")
    if publication_results.empty:
        st.info("No PubMed papers found for the selected metabolite-topic combinations.")
    else:
        st.dataframe(publication_results, use_container_width=True)

        st.download_button(
            "Download PubMed results CSV",
            data=publication_results.to_csv(index=False).encode("utf-8"),
            file_name="metabolite_pubmed_results.csv",
            mime="text/csv"
        )

    st.subheader("3. Paper-level summary")
    for metabolite in metabolites:
        st.markdown(f"### {metabolite}")
        sub = publication_results[publication_results.get("metabolite", "") == metabolite] if not publication_results.empty else pd.DataFrame()
        if sub.empty:
            st.write("No papers found.")
        else:
            for _, row in sub.head(5).iterrows():
                if "error" in row and pd.notna(row["error"]):
                    st.error(row["error"])
                    continue
                st.markdown(f"**{row.get('Title', '')}**")
                st.write(f"{row.get('Authors', '')} | {row.get('Journal', '')} | {row.get('Year', '')}")
                if row.get("PubMed_URL", ""):
                    st.markdown(f"[Open PubMed record]({row.get('PubMed_URL', '')})")
                abstract = row.get("Abstract", "")
                if abstract:
                    st.write(abstract[:1200] + ("..." if len(abstract) > 1200 else ""))
                st.divider()

# ----------------------------
# Database preview and template
# ----------------------------
with st.expander("Preview current annotation database"):
    st.dataframe(annotation_db, use_container_width=True)

with st.expander("Required annotation database format"):
    st.markdown("""
Your CSV should ideally contain these columns:

- metabolite
- synonyms
- compound_class
- pathway
- superpathway
- physiological_compartment
- hmdb_id
- kegg_id
- notes

You can start with the demo table, then expand it manually or from HMDB/KEGG/LipidMaps/ClassyFire-derived annotation exports.
""")

    template = pd.DataFrame(columns=[
        "metabolite", "synonyms", "compound_class", "pathway", "superpathway",
        "physiological_compartment", "hmdb_id", "kegg_id", "notes"
    ])

    st.download_button(
        "Download blank annotation database template",
        data=template.to_csv(index=False).encode("utf-8"),
        file_name="metabolite_annotation_template.csv",
        mime="text/csv"
    )
