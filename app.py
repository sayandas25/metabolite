
import streamlit as st
import pandas as pd
import requests
import xml.etree.ElementTree as ET
import re
import time
from urllib.parse import quote_plus

st.set_page_config(
    page_title="Metabolite Online Annotation & Literature Search",
    page_icon="🧬",
    layout="wide"
)

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
KEGG_BASE = "https://rest.kegg.jp"
NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

HEADERS = {
    "User-Agent": "MetaboliteOnlineAnnotationApp/1.0"
}


# =========================================================
# General helpers
# =========================================================
def safe_get(url, params=None, timeout=30):
    try:
        r = requests.get(url, params=params, headers=HEADERS, timeout=timeout)
        if r.status_code == 404:
            return None, "Not found"
        r.raise_for_status()
        return r, None
    except Exception as e:
        return None, str(e)


def split_metabolites(text):
    items = re.split(r"[\n,;]+", text)
    items = [x.strip() for x in items if x.strip()]
    return list(dict.fromkeys(items))


def first_nonempty(*vals):
    for v in vals:
        if v is not None and str(v).strip() != "":
            return str(v)
    return ""


def parse_kegg_flatfile(text):
    result = {}
    current = None
    for line in text.splitlines():
        if not line.strip():
            continue
        key = line[:12].strip()
        val = line[12:].strip()
        if key:
            current = key
            result[current] = val
        else:
            if current:
                result[current] += " " + val
    return result


# =========================================================
# Local database
# =========================================================
@st.cache_data(show_spinner=False)
def load_local_db(uploaded_file=None):
    if uploaded_file is not None:
        df = pd.read_csv(uploaded_file)
    else:
        df = pd.DataFrame([
            {
                "metabolite": "L-Cysteine",
                "synonyms": "Cysteine; Cys; L Cysteine",
                "compound_class": "Sulfur-containing amino acid",
                "pathway": "Cysteine and methionine metabolism; glutathione metabolism; taurine and hypotaurine metabolism",
                "superpathway": "Amino acid metabolism; redox metabolism",
                "physiological_compartment": "Plasma; cytosol; liver; placenta; immune cells",
                "hmdb_id": "HMDB0000574",
                "kegg_id": "C00097",
                "pubchem_cid": "5862",
                "inchikey": "XUJNEKJLAYXESH-REOHCLBHSA-N",
                "notes": "Important sulfur amino acid involved in redox balance, glutathione synthesis, and oxidative stress biology."
            },
            {
                "metabolite": "L-Arginine",
                "synonyms": "Arginine; Arg; L Arg",
                "compound_class": "Basic amino acid",
                "pathway": "Arginine and proline metabolism; nitric oxide biosynthesis; urea cycle",
                "superpathway": "Amino acid metabolism; vascular regulation",
                "physiological_compartment": "Plasma; cytosol; placenta; endothelial cells",
                "hmdb_id": "HMDB0000517",
                "kegg_id": "C00062",
                "pubchem_cid": "6322",
                "inchikey": "ODKSFYDXXFIFQN-BYPYZUCNSA-N",
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
                "pubchem_cid": "10917",
                "inchikey": "PHIQHXFUZVPYII-ZCFIWIBFSA-N",
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
                "pubchem_cid": "6305",
                "inchikey": "QIVBCDIJIAJPQS-VIFPVBQESA-N",
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
                "pubchem_cid": "305",
                "inchikey": "OEYIOHPDSNJKLS-UHFFFAOYSA-N",
                "notes": "Important for membrane phospholipids and methylation biology."
            }
        ])

    required_cols = [
        "metabolite", "synonyms", "compound_class", "pathway", "superpathway",
        "physiological_compartment", "hmdb_id", "kegg_id", "pubchem_cid",
        "inchikey", "notes"
    ]

    for col in required_cols:
        if col not in df.columns:
            df[col] = ""

    return df


def find_local_annotation(name, db):
    q = str(name).strip().lower()

    temp = db.copy()
    temp["metabolite_norm"] = temp["metabolite"].astype(str).str.lower().str.strip()
    temp["synonyms_norm"] = temp["synonyms"].astype(str).str.lower()

    exact = temp[temp["metabolite_norm"] == q]
    if not exact.empty:
        return exact.iloc[0].to_dict(), "Exact local match"

    syn = temp[temp["synonyms_norm"].str.contains(q, regex=False, na=False)]
    if not syn.empty:
        return syn.iloc[0].to_dict(), "Local synonym match"

    partial = temp[
        temp["metabolite_norm"].str.contains(q, regex=False, na=False) |
        temp["synonyms_norm"].str.contains(q, regex=False, na=False)
    ]
    if not partial.empty:
        return partial.iloc[0].to_dict(), "Partial local match"

    return {}, "No local match"


# =========================================================
# PubChem
# =========================================================
@st.cache_data(show_spinner=False, ttl=86400)
def pubchem_lookup(name):
    cid_url = f"{PUBCHEM_BASE}/compound/name/{quote_plus(name)}/cids/JSON"
    r, err = safe_get(cid_url)

    if r is None:
        return {
            "pubchem_status": "not_found",
            "pubchem_error": err
        }

    try:
        cids = r.json().get("IdentifierList", {}).get("CID", [])
        if not cids:
            return {"pubchem_status": "not_found", "pubchem_error": "No CID returned"}
        cid = str(cids[0])
    except Exception as e:
        return {"pubchem_status": "error", "pubchem_error": str(e)}

    prop_url = f"{PUBCHEM_BASE}/compound/cid/{cid}/property/IUPACName,MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES,InChIKey/JSON"
    pr, perr = safe_get(prop_url)

    props = {}
    if pr is not None:
        try:
            props = pr.json().get("PropertyTable", {}).get("Properties", [{}])[0]
        except Exception:
            props = {}

    syn_url = f"{PUBCHEM_BASE}/compound/cid/{cid}/synonyms/JSON"
    sr, serr = safe_get(syn_url)
    synonyms = []
    if sr is not None:
        try:
            synonyms = sr.json().get("InformationList", {}).get("Information", [{}])[0].get("Synonym", [])[:20]
        except Exception:
            synonyms = []

    return {
        "pubchem_status": "found",
        "pubchem_cid_online": cid,
        "pubchem_iupac_name": props.get("IUPACName", ""),
        "pubchem_formula": props.get("MolecularFormula", ""),
        "pubchem_molecular_weight": props.get("MolecularWeight", ""),
        "pubchem_canonical_smiles": props.get("CanonicalSMILES", ""),
        "pubchem_isomeric_smiles": props.get("IsomericSMILES", ""),
        "pubchem_inchikey": props.get("InChIKey", ""),
        "pubchem_synonyms": "; ".join(synonyms),
        "pubchem_url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}",
        "pubchem_error": ""
    }


# =========================================================
# KEGG
# =========================================================
@st.cache_data(show_spinner=False, ttl=86400)
def kegg_lookup_by_name(name):
    find_url = f"{KEGG_BASE}/find/compound/{quote_plus(name)}"
    r, err = safe_get(find_url)

    if r is None or not r.text.strip():
        return {
            "kegg_status": "not_found",
            "kegg_error": err or "No KEGG compound found"
        }

    first_line = r.text.strip().splitlines()[0]
    try:
        entry, label = first_line.split("\t", 1)
        kegg_id = entry.replace("cpd:", "")
    except Exception:
        return {
            "kegg_status": "error",
            "kegg_error": f"Could not parse KEGG result: {first_line}"
        }

    return kegg_lookup_by_id(kegg_id)


@st.cache_data(show_spinner=False, ttl=86400)
def kegg_lookup_by_id(kegg_id):
    if not kegg_id:
        return {"kegg_status": "not_found", "kegg_error": "No KEGG ID available"}

    kegg_id = str(kegg_id).replace("cpd:", "").strip()
    get_url = f"{KEGG_BASE}/get/cpd:{kegg_id}"
    r, err = safe_get(get_url)

    if r is None:
        return {
            "kegg_status": "not_found",
            "kegg_error": err,
            "kegg_id_online": kegg_id
        }

    parsed = parse_kegg_flatfile(r.text)

    link_url = f"{KEGG_BASE}/link/pathway/cpd:{kegg_id}"
    lr, lerr = safe_get(link_url)

    pathway_ids = []
    pathway_names = []

    if lr is not None and lr.text.strip():
        for line in lr.text.strip().splitlines():
            parts = line.split("\t")
            if len(parts) == 2:
                pathway_ids.append(parts[1].replace("path:", ""))

    for pid in pathway_ids[:10]:
        pr, _ = safe_get(f"{KEGG_BASE}/get/{pid}")
        if pr is not None:
            pparsed = parse_kegg_flatfile(pr.text)
            pname = pparsed.get("NAME", "")
            if pname:
                pathway_names.append(f"{pid}: {pname}")

    return {
        "kegg_status": "found",
        "kegg_id_online": kegg_id,
        "kegg_name": parsed.get("NAME", ""),
        "kegg_formula": parsed.get("FORMULA", ""),
        "kegg_exact_mass": parsed.get("EXACT_MASS", ""),
        "kegg_reaction": parsed.get("REACTION", ""),
        "kegg_enzyme": parsed.get("ENZYME", ""),
        "kegg_pathway_ids": "; ".join(pathway_ids),
        "kegg_pathway_names": "; ".join(pathway_names),
        "kegg_url": f"https://www.kegg.jp/entry/cpd:{kegg_id}",
        "kegg_error": ""
    }


# =========================================================
# ClassyFire
# =========================================================
@st.cache_data(show_spinner=False, ttl=86400)
def classyfire_lookup(inchikey):
    if not inchikey:
        return {
            "classyfire_status": "not_attempted",
            "classyfire_error": "No InChIKey available"
        }

    urls = [
        f"https://classyfire.wishartlab.com/entities/{inchikey}.json",
        f"http://classyfire.wishartlab.com/entities/{inchikey}.json"
    ]

    last_error = ""

    for url in urls:
        r, err = safe_get(url, timeout=30)
        if r is None:
            last_error = err
            continue

        try:
            data = r.json()

            def name_of(obj):
                if isinstance(obj, dict):
                    return obj.get("name", "")
                return ""

            return {
                "classyfire_status": "found",
                "classyfire_kingdom": name_of(data.get("kingdom")),
                "classyfire_superclass": name_of(data.get("superclass")),
                "classyfire_class": name_of(data.get("class")),
                "classyfire_subclass": name_of(data.get("subclass")),
                "classyfire_direct_parent": name_of(data.get("direct_parent")),
                "classyfire_error": ""
            }
        except Exception as e:
            last_error = str(e)

    return {
        "classyfire_status": "error",
        "classyfire_error": last_error
    }


# =========================================================
# PubMed and Europe PMC
# =========================================================
@st.cache_data(show_spinner=False, ttl=43200)
def pubmed_search(metabolite, relation_topic, max_results=10, email="", api_key=""):
    query = f'("{metabolite}"[Title/Abstract] OR "{metabolite}"[MeSH Terms]) AND ("{relation_topic}"[Title/Abstract] OR "{relation_topic}"[MeSH Terms])'

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

    r, err = safe_get(f"{NCBI_EUTILS}/esearch.fcgi", params=params, timeout=30)
    if r is None:
        return [{"metabolite": metabolite, "source": "PubMed", "error": err}]

    try:
        pmids = r.json().get("esearchresult", {}).get("idlist", [])
    except Exception as e:
        return [{"metabolite": metabolite, "source": "PubMed", "error": str(e)}]

    if not pmids:
        return []

    fr, ferr = safe_get(
        f"{NCBI_EUTILS}/efetch.fcgi",
        params={"db": "pubmed", "id": ",".join(pmids), "retmode": "xml"},
        timeout=40
    )

    if fr is None:
        return [{"metabolite": metabolite, "source": "PubMed", "error": ferr}]

    try:
        root = ET.fromstring(fr.content)
    except Exception as e:
        return [{"metabolite": metabolite, "source": "PubMed", "error": str(e)}]

    rows = []
    for article in root.findall(".//PubmedArticle"):
        pmid = article.findtext(".//PMID", default="")
        title_el = article.find(".//ArticleTitle")
        title = "".join(title_el.itertext()) if title_el is not None else ""
        journal = article.findtext(".//Journal/Title", default="")
        year = article.findtext(".//PubDate/Year", default="")
        if not year:
            year = article.findtext(".//PubDate/MedlineDate", default="")

        abstracts = []
        for ab in article.findall(".//Abstract/AbstractText"):
            abstracts.append("".join(ab.itertext()))

        authors = []
        for author in article.findall(".//Author")[:8]:
            last = author.findtext("LastName", default="")
            fore = author.findtext("ForeName", default="")
            nm = f"{fore} {last}".strip()
            if nm:
                authors.append(nm)

        rows.append({
            "metabolite": metabolite,
            "relation_topic": relation_topic,
            "source": "PubMed",
            "PMID": pmid,
            "Title": title,
            "Authors": "; ".join(authors),
            "Journal": journal,
            "Year": year,
            "Abstract": " ".join(abstracts),
            "URL": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else ""
        })

    return rows


@st.cache_data(show_spinner=False, ttl=43200)
def europe_pmc_search(metabolite, relation_topic, max_results=10):
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    params = {
        "query": f'"{metabolite}" AND "{relation_topic}"',
        "format": "json",
        "pageSize": max_results,
        "sort": "RELEVANCE"
    }

    r, err = safe_get(url, params=params, timeout=30)

    if r is None:
        return [{"metabolite": metabolite, "source": "Europe PMC", "error": err}]

    try:
        results = r.json().get("resultList", {}).get("result", [])
    except Exception as e:
        return [{"metabolite": metabolite, "source": "Europe PMC", "error": str(e)}]

    rows = []
    for x in results:
        rows.append({
            "metabolite": metabolite,
            "relation_topic": relation_topic,
            "source": "Europe PMC",
            "PMID": x.get("pmid", ""),
            "Title": x.get("title", ""),
            "Authors": x.get("authorString", ""),
            "Journal": x.get("journalTitle", ""),
            "Year": x.get("pubYear", ""),
            "DOI": x.get("doi", ""),
            "Cited_By_Count": x.get("citedByCount", ""),
            "Open_Access": x.get("isOpenAccess", ""),
            "Abstract": x.get("abstractText", ""),
            "URL": f"https://europepmc.org/article/MED/{x.get('pmid')}" if x.get("pmid") else ""
        })

    return rows


# =========================================================
# Integrated annotation
# =========================================================
def annotate_metabolite(name, local_db, use_pubchem, use_kegg, use_classyfire):
    local, local_status = find_local_annotation(name, local_db)

    pubchem = pubchem_lookup(name) if use_pubchem else {}
    local_kegg_id = str(local.get("kegg_id", "") or "").strip()

    if use_kegg:
        if local_kegg_id:
            kegg = kegg_lookup_by_id(local_kegg_id)
        else:
            kegg = kegg_lookup_by_name(name)
    else:
        kegg = {}

    inchikey = first_nonempty(
        local.get("inchikey", ""),
        pubchem.get("pubchem_inchikey", "")
    )

    classy = classyfire_lookup(inchikey) if use_classyfire else {}

    best_compound_class = first_nonempty(
        local.get("compound_class", ""),
        classy.get("classyfire_class", ""),
        classy.get("classyfire_superclass", ""),
        classy.get("classyfire_direct_parent", "")
    )

    best_pathway = first_nonempty(
        local.get("pathway", ""),
        kegg.get("kegg_pathway_names", "")
    )

    best_superpathway = first_nonempty(
        local.get("superpathway", ""),
        classy.get("classyfire_superclass", "")
    )

    best_compartment = first_nonempty(
        local.get("physiological_compartment", "")
    )

    annotation_note = ""
    if local_status == "No local match":
        annotation_note = "No local curated annotation. Online sources were searched. Physiological compartment usually requires HMDB/local curation."
    else:
        annotation_note = "Local curated annotation found. Online sources were also searched."

    row = {
        "input_metabolite": name,

        "best_compound_class": best_compound_class,
        "best_pathway": best_pathway,
        "best_superpathway": best_superpathway,
        "best_physiological_compartment": best_compartment,
        "annotation_note": annotation_note,

        "local_status": local_status,
        "local_matched_metabolite": local.get("metabolite", ""),
        "hmdb_id": local.get("hmdb_id", ""),
        "local_kegg_id": local.get("kegg_id", ""),
        "local_pubchem_cid": local.get("pubchem_cid", ""),

        "pubchem_status": pubchem.get("pubchem_status", ""),
        "pubchem_cid": first_nonempty(local.get("pubchem_cid", ""), pubchem.get("pubchem_cid_online", "")),
        "pubchem_formula": pubchem.get("pubchem_formula", ""),
        "pubchem_molecular_weight": pubchem.get("pubchem_molecular_weight", ""),
        "pubchem_inchikey": pubchem.get("pubchem_inchikey", ""),
        "pubchem_iupac_name": pubchem.get("pubchem_iupac_name", ""),
        "pubchem_url": pubchem.get("pubchem_url", ""),
        "pubchem_error": pubchem.get("pubchem_error", ""),

        "kegg_status": kegg.get("kegg_status", ""),
        "kegg_id": first_nonempty(local.get("kegg_id", ""), kegg.get("kegg_id_online", "")),
        "kegg_name": kegg.get("kegg_name", ""),
        "kegg_formula": kegg.get("kegg_formula", ""),
        "kegg_exact_mass": kegg.get("kegg_exact_mass", ""),
        "kegg_pathway_ids": kegg.get("kegg_pathway_ids", ""),
        "kegg_pathway_names": kegg.get("kegg_pathway_names", ""),
        "kegg_url": kegg.get("kegg_url", ""),
        "kegg_error": kegg.get("kegg_error", ""),

        "classyfire_status": classy.get("classyfire_status", ""),
        "classyfire_kingdom": classy.get("classyfire_kingdom", ""),
        "classyfire_superclass": classy.get("classyfire_superclass", ""),
        "classyfire_class": classy.get("classyfire_class", ""),
        "classyfire_subclass": classy.get("classyfire_subclass", ""),
        "classyfire_direct_parent": classy.get("classyfire_direct_parent", ""),
        "classyfire_error": classy.get("classyfire_error", "")
    }

    return row


# =========================================================
# UI
# =========================================================
st.title("🧬 Metabolite Online Annotation & Literature Search")

st.markdown(
    """
Enter one metabolite or a list of metabolites. The app searches a local curated table plus online sources:
**PubChem**, **KEGG**, **ClassyFire**, **PubMed**, and **Europe PMC**.
"""
)

with st.sidebar:
    st.header("Settings")

    uploaded_db = st.file_uploader(
        "Optional: upload curated annotation CSV",
        type=["csv"]
    )

    local_db = load_local_db(uploaded_db)

    st.markdown("### Annotation sources")
    use_pubchem = st.checkbox("Search PubChem", value=True)
    use_kegg = st.checkbox("Search KEGG", value=True)
    use_classyfire = st.checkbox("Search ClassyFire", value=True)

    st.markdown("### Literature sources")
    use_pubmed = st.checkbox("Search PubMed", value=True)
    use_epmc = st.checkbox("Search Europe PMC", value=True)

    max_papers = st.slider("Maximum papers per metabolite/source", 3, 30, 10)

    ncbi_email = st.text_input("NCBI email, optional", placeholder="your.email@example.com")
    ncbi_api_key = st.text_input("NCBI API key, optional", type="password")

tab1, tab2, tab3 = st.tabs(["Run app", "Local database", "Template"])

with tab1:
    col1, col2 = st.columns([1, 1])

    with col1:
        relation_topic = st.text_input(
            "Relation topic for literature search",
            placeholder="Example: preeclampsia, gestational hypertension, preterm birth"
        )

    with col2:
        input_mode = st.radio(
            "Input mode",
            ["Single metabolite", "List of metabolites"],
            horizontal=True
        )

    if input_mode == "Single metabolite":
        metabolite_text = st.text_input("Metabolite name", placeholder="Example: Cysteine")
    else:
        metabolite_text = st.text_area("Metabolite list", placeholder="Cysteine\nArginine\nCarnitine")

    run = st.button("Run online annotation and literature search", type="primary")

    if run:
        if not metabolite_text.strip():
            st.error("Please enter at least one metabolite.")
            st.stop()

        metabolites = split_metabolites(metabolite_text)

        annotation_rows = []
        literature_rows = []

        progress = st.progress(0)
        status = st.empty()

        for i, met in enumerate(metabolites):
            status.write(f"Searching online databases for: {met}")

            ann = annotate_metabolite(
                met,
                local_db=local_db,
                use_pubchem=use_pubchem,
                use_kegg=use_kegg,
                use_classyfire=use_classyfire
            )
            annotation_rows.append(ann)

            if relation_topic.strip():
                if use_pubmed:
                    literature_rows.extend(pubmed_search(
                        met,
                        relation_topic,
                        max_results=max_papers,
                        email=ncbi_email,
                        api_key=ncbi_api_key
                    ))
                    time.sleep(0.34)

                if use_epmc:
                    literature_rows.extend(europe_pmc_search(
                        met,
                        relation_topic,
                        max_results=max_papers
                    ))
                    time.sleep(0.34)

            progress.progress((i + 1) / len(metabolites))

        status.write("Completed.")

        ann_df = pd.DataFrame(annotation_rows)
        lit_df = pd.DataFrame(literature_rows)

        st.subheader("1. Integrated annotation results")
        st.dataframe(ann_df, use_container_width=True)

        st.download_button(
            "Download annotation CSV",
            data=ann_df.to_csv(index=False).encode("utf-8"),
            file_name="integrated_metabolite_annotation.csv",
            mime="text/csv"
        )

        if relation_topic.strip():
            st.subheader("2. Literature results")
            if lit_df.empty:
                st.info("No literature results found for this metabolite-topic combination.")
            else:
                st.dataframe(lit_df, use_container_width=True)

                st.download_button(
                    "Download literature CSV",
                    data=lit_df.to_csv(index=False).encode("utf-8"),
                    file_name="metabolite_literature_results.csv",
                    mime="text/csv"
                )

                st.subheader("3. Paper cards")
                for _, row in lit_df.head(20).iterrows():
                    if pd.notna(row.get("error", None)):
                        st.error(row.get("error"))
                        continue

                    st.markdown(f"**{row.get('Title', '')}**")
                    st.write(f"{row.get('source', '')} | {row.get('Journal', '')} | {row.get('Year', '')}")
                    if row.get("URL", ""):
                        st.markdown(f"[Open record]({row.get('URL')})")
                    if row.get("Abstract", ""):
                        abstract = str(row.get("Abstract", ""))
                        st.write(abstract[:1200] + ("..." if len(abstract) > 1200 else ""))
                    st.divider()
        else:
            st.info("No relation topic was entered, so literature search was skipped.")

with tab2:
    st.subheader("Current local curated database")
    st.dataframe(local_db, use_container_width=True)

with tab3:
    st.subheader("Annotation CSV template")
    template = pd.DataFrame(columns=[
        "metabolite",
        "synonyms",
        "compound_class",
        "pathway",
        "superpathway",
        "physiological_compartment",
        "hmdb_id",
        "kegg_id",
        "pubchem_cid",
        "inchikey",
        "notes"
    ])
    st.dataframe(template, use_container_width=True)

    st.download_button(
        "Download template CSV",
        data=template.to_csv(index=False).encode("utf-8"),
        file_name="metabolite_annotation_template.csv",
        mime="text/csv"
    )
