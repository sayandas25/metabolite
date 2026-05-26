
import streamlit as st
import pandas as pd
import requests
import xml.etree.ElementTree as ET
import re
import time
from urllib.parse import quote_plus

st.set_page_config(
    page_title="Metabolite Annotation v4 - No Blank Outputs",
    page_icon="🧬",
    layout="wide"
)

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
KEGG_BASE = "https://rest.kegg.jp"
NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OLS_BASE = "https://www.ebi.ac.uk/ols/api"

HEADERS = {"User-Agent": "MetaboliteAnnotationNoBlankV4/1.0"}

# -----------------------------
# Basic helpers
# -----------------------------
def safe_get(url, params=None, timeout=30):
    try:
        r = requests.get(url, params=params, headers=HEADERS, timeout=timeout)
        if r.status_code == 404:
            return None, "Not found"
        r.raise_for_status()
        return r, None
    except Exception as e:
        return None, str(e)

def clean(x):
    if x is None:
        return ""
    return str(x).strip()

def is_blank(x):
    return x is None or str(x).strip() in ["", "nan", "None", "not available", "Not available"]

def first_value(*vals, default="Not available"):
    for v in vals:
        if not is_blank(v):
            return str(v).strip()
    return default

def split_metabolites(text):
    items = re.split(r"[\n,;]+", text)
    return list(dict.fromkeys([x.strip() for x in items if x.strip()]))

def name_variants(name):
    n = clean(name)
    variants = []

    def add(x):
        x = clean(x)
        if x and x not in variants:
            variants.append(x)

    add(n)
    add(n.replace("_", " "))
    add(n.replace("-", " "))
    add(re.sub(r"\s+", " ", n))
    add(n.replace("(", "").replace(")", ""))
    add(n.replace("[", "").replace("]", ""))

    low = n.lower()
    if not low.startswith("l-"):
        add("L-" + n)
        add("L " + n)
    if low.startswith("l-"):
        add(n[2:])
    if low.startswith("l "):
        add(n[2:])

    add(n.replace("α", "alpha").replace("β", "beta").replace("γ", "gamma"))
    add(n.replace("alpha", "α").replace("beta", "β").replace("gamma", "γ"))

    # Remove common LC-MS annotations
    cleaned = re.sub(r"\bpos\b|\bneg\b|\bpositive\b|\bnegative\b", "", n, flags=re.I)
    cleaned = re.sub(r"\[.*?\]", "", cleaned)
    cleaned = re.sub(r"\s+", " ", cleaned).strip()
    add(cleaned)

    return variants

def parse_kegg(text):
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
        elif current:
            result[current] += " " + val
    return result

# -----------------------------
# Local curated database
# -----------------------------
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
                "physiological_compartment": "Blood/plasma; cytosol; liver; placenta; immune cells",
                "hmdb_id": "HMDB0000574",
                "kegg_id": "C00097",
                "pubchem_cid": "5862",
                "inchikey": "XUJNEKJLAYXESH-REOHCLBHSA-N",
                "notes": "Sulfur amino acid involved in glutathione synthesis and redox balance."
            },
            {
                "metabolite": "L-Arginine",
                "synonyms": "Arginine; Arg; L Arg",
                "compound_class": "Basic amino acid",
                "pathway": "Arginine and proline metabolism; nitric oxide biosynthesis; urea cycle",
                "superpathway": "Amino acid metabolism; vascular regulation",
                "physiological_compartment": "Blood/plasma; cytosol; placenta; endothelial cells; liver",
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
                "physiological_compartment": "Blood/plasma; mitochondria; muscle; liver; placenta",
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
                "physiological_compartment": "Blood/plasma; gut; placenta; immune cells; brain",
                "hmdb_id": "HMDB0000929",
                "kegg_id": "C00078",
                "pubchem_cid": "6305",
                "inchikey": "QIVBCDIJIAJPQS-VIFPVBQESA-N",
                "notes": "Precursor of kynurenine, serotonin, melatonin and microbial indoles."
            },
            {
                "metabolite": "Choline",
                "synonyms": "2-Hydroxyethyltrimethylammonium; Bilineurine",
                "compound_class": "Quaternary ammonium compound",
                "pathway": "Phosphatidylcholine metabolism; one-carbon metabolism",
                "superpathway": "Lipid metabolism; methyl donor metabolism",
                "physiological_compartment": "Blood/plasma; liver; placenta; cell membrane; brain",
                "hmdb_id": "HMDB0000097",
                "kegg_id": "C00114",
                "pubchem_cid": "305",
                "inchikey": "OEYIOHPDSNJKLS-UHFFFAOYSA-N",
                "notes": "Important for membrane phospholipids and methylation biology."
            }
        ])

    required = [
        "metabolite", "synonyms", "compound_class", "pathway", "superpathway",
        "physiological_compartment", "hmdb_id", "kegg_id", "pubchem_cid",
        "inchikey", "notes"
    ]
    for c in required:
        if c not in df.columns:
            df[c] = ""
    return df

def find_local(name, db):
    variants = [v.lower() for v in name_variants(name)]
    temp = db.copy()
    temp["metabolite_norm"] = temp["metabolite"].astype(str).str.lower().str.strip()
    temp["synonyms_norm"] = temp["synonyms"].astype(str).str.lower()

    for q in variants:
        hit = temp[temp["metabolite_norm"] == q]
        if not hit.empty:
            return hit.iloc[0].to_dict(), "Exact local match"

    for q in variants:
        hit = temp[temp["synonyms_norm"].str.contains(q, regex=False, na=False)]
        if not hit.empty:
            return hit.iloc[0].to_dict(), "Local synonym match"

    return {}, "No local match"

# -----------------------------
# Online databases
# -----------------------------
@st.cache_data(show_spinner=False, ttl=86400)
def pubchem_lookup(name):
    errors = []
    for v in name_variants(name):
        r, err = safe_get(f"{PUBCHEM_BASE}/compound/name/{quote_plus(v)}/cids/JSON")
        if r is None:
            errors.append(f"{v}: {err}")
            continue
        try:
            cids = r.json().get("IdentifierList", {}).get("CID", [])
            if not cids:
                errors.append(f"{v}: no CID")
                continue
            cid = str(cids[0])
        except Exception as e:
            errors.append(f"{v}: {e}")
            continue

        props = {}
        pr, _ = safe_get(
            f"{PUBCHEM_BASE}/compound/cid/{cid}/property/"
            "IUPACName,MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES,InChIKey/JSON"
        )
        if pr is not None:
            try:
                props = pr.json().get("PropertyTable", {}).get("Properties", [{}])[0]
            except Exception:
                props = {}

        synonyms = []
        sr, _ = safe_get(f"{PUBCHEM_BASE}/compound/cid/{cid}/synonyms/JSON")
        if sr is not None:
            try:
                synonyms = sr.json().get("InformationList", {}).get("Information", [{}])[0].get("Synonym", [])[:25]
            except Exception:
                synonyms = []

        return {
            "pubchem_status": "found",
            "pubchem_lookup_name": v,
            "pubchem_cid": cid,
            "pubchem_iupac_name": props.get("IUPACName", ""),
            "pubchem_formula": props.get("MolecularFormula", ""),
            "pubchem_molecular_weight": props.get("MolecularWeight", ""),
            "pubchem_smiles": props.get("IsomericSMILES", "") or props.get("CanonicalSMILES", ""),
            "pubchem_inchikey": props.get("InChIKey", ""),
            "pubchem_synonyms": "; ".join(synonyms),
            "pubchem_url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}",
            "pubchem_error": ""
        }

    return {
        "pubchem_status": "not_found",
        "pubchem_lookup_name": "Tried: " + "; ".join(name_variants(name)[:5]),
        "pubchem_cid": "",
        "pubchem_iupac_name": "",
        "pubchem_formula": "",
        "pubchem_molecular_weight": "",
        "pubchem_smiles": "",
        "pubchem_inchikey": "",
        "pubchem_synonyms": "",
        "pubchem_url": "",
        "pubchem_error": " | ".join(errors[:3]) if errors else "No PubChem result"
    }

@st.cache_data(show_spinner=False, ttl=86400)
def kegg_lookup_by_id(kegg_id):
    kegg_id = clean(kegg_id).replace("cpd:", "")
    if not kegg_id:
        return {"kegg_status": "not_found"}

    r, err = safe_get(f"{KEGG_BASE}/get/cpd:{kegg_id}")
    if r is None:
        return {"kegg_status": "not_found", "kegg_error": err, "kegg_id": kegg_id}

    parsed = parse_kegg(r.text)

    pathway_ids = []
    pathway_names = []
    lr, _ = safe_get(f"{KEGG_BASE}/link/pathway/cpd:{kegg_id}")
    if lr is not None and lr.text.strip():
        for line in lr.text.strip().splitlines():
            parts = line.split("\t")
            if len(parts) == 2:
                pathway_ids.append(parts[1].replace("path:", ""))

    for pid in pathway_ids[:12]:
        pr, _ = safe_get(f"{KEGG_BASE}/get/{pid}")
        if pr is not None:
            pp = parse_kegg(pr.text)
            if pp.get("NAME"):
                pathway_names.append(f"{pid}: {pp.get('NAME')}")

    return {
        "kegg_status": "found",
        "kegg_id": kegg_id,
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

@st.cache_data(show_spinner=False, ttl=86400)
def kegg_lookup_by_name(name):
    errors = []
    for v in name_variants(name):
        r, err = safe_get(f"{KEGG_BASE}/find/compound/{quote_plus(v)}")
        if r is None or not r.text.strip():
            errors.append(f"{v}: {err or 'no result'}")
            continue
        line = r.text.strip().splitlines()[0]
        try:
            entry = line.split("\t", 1)[0]
            kid = entry.replace("cpd:", "")
            out = kegg_lookup_by_id(kid)
            out["kegg_lookup_name"] = v
            return out
        except Exception:
            errors.append(f"{v}: parse failed")
    return {
        "kegg_status": "not_found",
        "kegg_lookup_name": "Tried: " + "; ".join(name_variants(name)[:5]),
        "kegg_id": "",
        "kegg_name": "",
        "kegg_formula": "",
        "kegg_exact_mass": "",
        "kegg_reaction": "",
        "kegg_enzyme": "",
        "kegg_pathway_ids": "",
        "kegg_pathway_names": "",
        "kegg_url": "",
        "kegg_error": " | ".join(errors[:3]) if errors else "No KEGG result"
    }

@st.cache_data(show_spinner=False, ttl=86400)
def chebi_lookup(name):
    errors = []
    for v in name_variants(name):
        for exact in ["true", "false"]:
            r, err = safe_get(
                f"{OLS_BASE}/search",
                params={"q": v, "ontology": "chebi", "rows": 5, "exact": exact},
                timeout=30
            )
            if r is None:
                errors.append(f"{v}: {err}")
                continue
            try:
                docs = r.json().get("response", {}).get("docs", [])
            except Exception as e:
                errors.append(f"{v}: {e}")
                continue
            if docs:
                d = docs[0]
                desc = d.get("description", "")
                if isinstance(desc, list):
                    desc = desc[0] if desc else ""
                return {
                    "chebi_status": "found",
                    "chebi_lookup_name": v,
                    "chebi_label": d.get("label", ""),
                    "chebi_short_form": d.get("short_form", ""),
                    "chebi_iri": d.get("iri", ""),
                    "chebi_description": desc,
                    "chebi_error": ""
                }

    return {
        "chebi_status": "not_found",
        "chebi_lookup_name": "Tried: " + "; ".join(name_variants(name)[:5]),
        "chebi_label": "",
        "chebi_short_form": "",
        "chebi_iri": "",
        "chebi_description": "",
        "chebi_error": " | ".join(errors[:3]) if errors else "No ChEBI result"
    }

@st.cache_data(show_spinner=False, ttl=86400)
def classyfire_lookup(inchikey):
    if is_blank(inchikey):
        return {
            "classyfire_status": "not_attempted",
            "classyfire_kingdom": "",
            "classyfire_superclass": "",
            "classyfire_class": "",
            "classyfire_subclass": "",
            "classyfire_direct_parent": "",
            "classyfire_error": "No InChIKey available"
        }

    last_err = ""
    for url in [
        f"https://classyfire.wishartlab.com/entities/{inchikey}.json",
        f"http://classyfire.wishartlab.com/entities/{inchikey}.json"
    ]:
        r, err = safe_get(url, timeout=30)
        if r is None:
            last_err = err
            continue
        try:
            data = r.json()
            def pick(obj):
                return obj.get("name", "") if isinstance(obj, dict) else ""
            return {
                "classyfire_status": "found",
                "classyfire_kingdom": pick(data.get("kingdom")),
                "classyfire_superclass": pick(data.get("superclass")),
                "classyfire_class": pick(data.get("class")),
                "classyfire_subclass": pick(data.get("subclass")),
                "classyfire_direct_parent": pick(data.get("direct_parent")),
                "classyfire_error": ""
            }
        except Exception as e:
            last_err = str(e)

    return {
        "classyfire_status": "error",
        "classyfire_kingdom": "",
        "classyfire_superclass": "",
        "classyfire_class": "",
        "classyfire_subclass": "",
        "classyfire_direct_parent": "",
        "classyfire_error": last_err
    }

# -----------------------------
# Rule-based fallback
# -----------------------------
def infer_fallback(name, pathway="", class_text="", synonyms="", formula=""):
    text = " ".join([name, pathway, class_text, synonyms, formula]).lower()

    cls = "Likely exogenous/xenobiotic or unclassified small molecule"
    superpath = "Exogenous/xenobiotic metabolism or unclassified metabolism"
    path = "Unclassified or database-specific pathway"
    compartment = "Likely blood/plasma or exposure-related compartment"
    confidence = "Low"

    amino = ["alanine","arginine","asparagine","aspartate","cysteine","glutamate","glutamine","glycine",
             "histidine","isoleucine","leucine","lysine","methionine","phenylalanine","proline","serine",
             "threonine","tryptophan","tyrosine","valine","ornithine","citrulline","taurine"]
    if any(a in text for a in amino) or "amino acid" in text:
        cls = "Amino acid or amino-acid derivative"
        superpath = "Amino acid metabolism"
        path = first_value(pathway, "Amino acid metabolism")
        compartment = "Blood/plasma; liver; cytosol; placenta/fetal-maternal interface if pregnancy-related"
        confidence = "Moderate"

    if any(x in text for x in ["cysteine", "glutathione", "methionine", "taurine", "sulfur"]):
        cls = "Sulfur-containing amino acid or redox metabolite"
        superpath = "Amino acid metabolism; redox metabolism"
        path = first_value(pathway, "Cysteine, methionine, glutathione or taurine metabolism")
        compartment = "Blood/plasma; liver; cytosol; placenta; immune cells"
        confidence = "Moderate"

    if any(x in text for x in ["carnitine", "acylcarnitine", "palmitoyl", "oleoyl"]):
        cls = "Carnitine or acylcarnitine lipid-energy metabolite"
        superpath = "Fatty acid oxidation; mitochondrial energy metabolism"
        path = first_value(pathway, "Carnitine shuttle and fatty acid beta-oxidation")
        compartment = "Blood/plasma; mitochondria; liver; skeletal muscle; placenta"
        confidence = "Moderate"

    if any(x in text for x in ["phosphatidyl", "lysophosph", "ceramide", "sphingo", "sphingomyelin", "lipid", "pc(", "pe("]):
        cls = "Lipid or membrane lipid species"
        superpath = "Lipid metabolism; membrane remodeling"
        path = first_value(pathway, "Glycerophospholipid or sphingolipid metabolism")
        compartment = "Blood/plasma lipoproteins; cell membrane; liver; placenta"
        confidence = "Moderate"

    if any(x in text for x in ["cholate","deoxycholate","chenodeoxycholate","taurocholate","glycocholate","bile acid"]):
        cls = "Bile acid or bile acid conjugate"
        superpath = "Steroid and bile acid metabolism"
        path = first_value(pathway, "Primary/secondary bile acid metabolism")
        compartment = "Liver; bile; intestine; blood/plasma"
        confidence = "Moderate"

    if any(x in text for x in ["steroid","estradiol","estrone","progesterone","cortisol","cholesterol","androgen"]):
        cls = "Steroid or steroid-derived metabolite"
        superpath = "Steroid hormone metabolism"
        path = first_value(pathway, "Steroid biosynthesis or steroid hormone metabolism")
        compartment = "Placenta; adrenal/gonadal tissue; liver; blood/plasma"
        confidence = "Moderate"

    if any(x in text for x in ["glucose","fructose","lactate","pyruvate","citrate","succinate","fumarate","malate"]):
        cls = "Carbohydrate or central carbon metabolism intermediate"
        superpath = "Carbohydrate and energy metabolism"
        path = first_value(pathway, "Glycolysis, TCA cycle or central carbon metabolism")
        compartment = "Blood/plasma; cytosol; mitochondria; liver; placenta"
        confidence = "Moderate"

    if any(x in text for x in ["uridine","adenosine","guanosine","inosine","hypoxanthine","xanthine","uric acid","purine","pyrimidine"]):
        cls = "Nucleotide, nucleoside or purine/pyrimidine derivative"
        superpath = "Nucleotide metabolism"
        path = first_value(pathway, "Purine or pyrimidine metabolism")
        compartment = "Blood/plasma; cytosol; kidney; liver; placenta"
        confidence = "Moderate"

    if any(x in text for x in ["indole","hippurate","cresol","butyrate","propionate","microbial"]):
        cls = "Microbial-host co-metabolite"
        superpath = "Gut microbiome-related metabolism"
        path = first_value(pathway, "Microbial aromatic amino acid or short-chain fatty acid metabolism")
        compartment = "Gut/intestine; blood/plasma; liver"
        confidence = "Moderate"

    if any(x in text for x in ["drug","xenobiotic","pesticide","phthalate","paraben","caffeine","nicotine","benzoate","glucuronide"]):
        cls = "Likely exogenous/xenobiotic compound or conjugated exposure marker"
        superpath = "Xenobiotic metabolism"
        path = first_value(pathway, "Phase I/II detoxification, sulfation or glucuronidation")
        compartment = "Exposure-related; liver; kidney; blood/plasma; urine"
        confidence = "Moderate"

    return cls, path, superpath, compartment, confidence

def no_blank_row(row):
    # first replace raw blanks
    for k in list(row.keys()):
        if is_blank(row[k]):
            row[k] = "Not available"

    cls_fb, path_fb, super_fb, comp_fb, conf = infer_fallback(
        row.get("input_metabolite", ""),
        pathway=first_value(row.get("pathway_local", ""), row.get("kegg_pathway_names", ""), default=""),
        class_text=first_value(row.get("compound_class_local", ""), row.get("classyfire_class", ""), row.get("classyfire_superclass", ""), default=""),
        synonyms=first_value(row.get("pubchem_synonyms", ""), default=""),
        formula=first_value(row.get("pubchem_formula", ""), default="")
    )

    row["best_compound_class"] = first_value(
        row.get("compound_class_local", ""),
        row.get("classyfire_class", ""),
        row.get("classyfire_superclass", ""),
        row.get("classyfire_direct_parent", ""),
        row.get("chebi_label", ""),
        cls_fb
    )
    row["best_pathway"] = first_value(row.get("pathway_local", ""), row.get("kegg_pathway_names", ""), path_fb)
    row["best_superpathway"] = first_value(row.get("superpathway_local", ""), row.get("classyfire_superclass", ""), super_fb)
    row["best_physiological_compartment"] = first_value(row.get("physiological_compartment_local", ""), comp_fb)
    row["compartment_source"] = "Local curated table" if row.get("physiological_compartment_local") != "Not available" else "Rule-based closest organ/compartment inference"
    row["fallback_confidence"] = conf

    for k in list(row.keys()):
        if is_blank(row[k]):
            row[k] = "Not available"
    return row

def annotate_one(name, db, use_pubchem=True, use_kegg=True, use_chebi=True, use_classyfire=True):
    local, local_status = find_local(name, db)
    pub = pubchem_lookup(name) if use_pubchem else {"pubchem_status": "not used"}
    kid = first_value(local.get("kegg_id", ""), default="")
    kegg = kegg_lookup_by_id(kid) if use_kegg and kid else (kegg_lookup_by_name(name) if use_kegg else {"kegg_status": "not used"})
    chebi = chebi_lookup(name) if use_chebi else {"chebi_status": "not used"}
    inchikey = first_value(local.get("inchikey", ""), pub.get("pubchem_inchikey", ""), default="")
    classy = classyfire_lookup(inchikey) if use_classyfire else {"classyfire_status": "not used"}

    row = {
        "input_metabolite": name,
        "local_status": local_status,
        "local_matched_metabolite": local.get("metabolite", ""),
        "compound_class_local": local.get("compound_class", ""),
        "pathway_local": local.get("pathway", ""),
        "superpathway_local": local.get("superpathway", ""),
        "physiological_compartment_local": local.get("physiological_compartment", ""),
        "hmdb_id": local.get("hmdb_id", ""),
        "local_kegg_id": local.get("kegg_id", ""),
        "local_pubchem_cid": local.get("pubchem_cid", ""),
        "local_notes": local.get("notes", ""),

        "pubchem_status": pub.get("pubchem_status", ""),
        "pubchem_lookup_name": pub.get("pubchem_lookup_name", ""),
        "pubchem_cid": first_value(local.get("pubchem_cid", ""), pub.get("pubchem_cid", ""), default=""),
        "pubchem_formula": pub.get("pubchem_formula", ""),
        "pubchem_molecular_weight": pub.get("pubchem_molecular_weight", ""),
        "pubchem_smiles": pub.get("pubchem_smiles", ""),
        "pubchem_inchikey": pub.get("pubchem_inchikey", ""),
        "pubchem_iupac_name": pub.get("pubchem_iupac_name", ""),
        "pubchem_synonyms": pub.get("pubchem_synonyms", ""),
        "pubchem_url": pub.get("pubchem_url", ""),
        "pubchem_error": pub.get("pubchem_error", ""),

        "kegg_status": kegg.get("kegg_status", ""),
        "kegg_id": first_value(local.get("kegg_id", ""), kegg.get("kegg_id", ""), default=""),
        "kegg_name": kegg.get("kegg_name", ""),
        "kegg_formula": kegg.get("kegg_formula", ""),
        "kegg_exact_mass": kegg.get("kegg_exact_mass", ""),
        "kegg_reaction": kegg.get("kegg_reaction", ""),
        "kegg_enzyme": kegg.get("kegg_enzyme", ""),
        "kegg_pathway_ids": kegg.get("kegg_pathway_ids", ""),
        "kegg_pathway_names": kegg.get("kegg_pathway_names", ""),
        "kegg_url": kegg.get("kegg_url", ""),
        "kegg_error": kegg.get("kegg_error", ""),

        "chebi_status": chebi.get("chebi_status", ""),
        "chebi_lookup_name": chebi.get("chebi_lookup_name", ""),
        "chebi_label": chebi.get("chebi_label", ""),
        "chebi_short_form": chebi.get("chebi_short_form", ""),
        "chebi_iri": chebi.get("chebi_iri", ""),
        "chebi_description": chebi.get("chebi_description", ""),
        "chebi_error": chebi.get("chebi_error", ""),

        "classyfire_status": classy.get("classyfire_status", ""),
        "classyfire_kingdom": classy.get("classyfire_kingdom", ""),
        "classyfire_superclass": classy.get("classyfire_superclass", ""),
        "classyfire_class": classy.get("classyfire_class", ""),
        "classyfire_subclass": classy.get("classyfire_subclass", ""),
        "classyfire_direct_parent": classy.get("classyfire_direct_parent", ""),
        "classyfire_error": classy.get("classyfire_error", "")
    }
    return no_blank_row(row)

# -----------------------------
# Literature search
# -----------------------------
@st.cache_data(show_spinner=False, ttl=43200)
def pubmed_search(metabolite, topic, max_results=10, email="", api_key=""):
    query = f'("{metabolite}"[Title/Abstract] OR "{metabolite}"[MeSH Terms]) AND ("{topic}"[Title/Abstract] OR "{topic}"[MeSH Terms])'
    params = {"db": "pubmed", "term": query, "retmode": "json", "retmax": max_results, "sort": "relevance"}
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key

    r, err = safe_get(f"{NCBI_EUTILS}/esearch.fcgi", params=params)
    if r is None:
        return [{"metabolite": metabolite, "source": "PubMed", "error": err}]
    try:
        pmids = r.json().get("esearchresult", {}).get("idlist", [])
    except Exception as e:
        return [{"metabolite": metabolite, "source": "PubMed", "error": str(e)}]
    if not pmids:
        return []

    fr, ferr = safe_get(f"{NCBI_EUTILS}/efetch.fcgi", params={"db": "pubmed", "id": ",".join(pmids), "retmode": "xml"})
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
        title = "".join(title_el.itertext()) if title_el is not None else "Not available"
        journal = article.findtext(".//Journal/Title", default="Not available")
        year = article.findtext(".//PubDate/Year", default="")
        if not year:
            year = article.findtext(".//PubDate/MedlineDate", default="Not available")

        abstracts = ["".join(ab.itertext()) for ab in article.findall(".//Abstract/AbstractText")]
        authors = []
        for a in article.findall(".//Author")[:8]:
            nm = f"{a.findtext('ForeName', default='')} {a.findtext('LastName', default='')}".strip()
            if nm:
                authors.append(nm)

        rows.append({
            "metabolite": metabolite,
            "relation_topic": topic,
            "source": "PubMed",
            "PMID": pmid or "Not available",
            "Title": title or "Not available",
            "Authors": "; ".join(authors) or "Not available",
            "Journal": journal or "Not available",
            "Year": year or "Not available",
            "Abstract": " ".join(abstracts) or "Not available",
            "URL": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "Not available"
        })
    return rows

@st.cache_data(show_spinner=False, ttl=43200)
def europe_pmc_search(metabolite, topic, max_results=10):
    r, err = safe_get(
        "https://www.ebi.ac.uk/europepmc/webservices/rest/search",
        params={"query": f'"{metabolite}" AND "{topic}"', "format": "json", "pageSize": max_results, "sort": "RELEVANCE"}
    )
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
            "relation_topic": topic,
            "source": "Europe PMC",
            "PMID": x.get("pmid", "Not available") or "Not available",
            "Title": x.get("title", "Not available") or "Not available",
            "Authors": x.get("authorString", "Not available") or "Not available",
            "Journal": x.get("journalTitle", "Not available") or "Not available",
            "Year": x.get("pubYear", "Not available") or "Not available",
            "DOI": x.get("doi", "Not available") or "Not available",
            "Cited_By_Count": x.get("citedByCount", "Not available") or "Not available",
            "Open_Access": x.get("isOpenAccess", "Not available") or "Not available",
            "Abstract": x.get("abstractText", "Not available") or "Not available",
            "URL": f"https://europepmc.org/article/MED/{x.get('pmid')}" if x.get("pmid") else "Not available"
        })
    return rows

# -----------------------------
# UI
# -----------------------------
st.title("🧬 Metabolite Annotation v4")
st.caption("No-blank output version: local database + PubChem + KEGG + ChEBI/OLS + ClassyFire + PubMed + Europe PMC")

with st.sidebar:
    st.header("Settings")
    uploaded_db = st.file_uploader("Optional curated annotation CSV", type=["csv"])
    local_db = load_local_db(uploaded_db)

    st.subheader("Annotation sources")
    use_pubchem = st.checkbox("Search PubChem", value=True)
    use_kegg = st.checkbox("Search KEGG", value=True)
    use_chebi = st.checkbox("Search ChEBI via EBI OLS", value=True)
    use_classyfire = st.checkbox("Search ClassyFire", value=True)

    st.subheader("Literature sources")
    use_pubmed = st.checkbox("Search PubMed", value=True)
    use_epmc = st.checkbox("Search Europe PMC", value=True)
    max_papers = st.slider("Maximum papers per source", 3, 30, 10)

    st.subheader("NCBI")
    ncbi_email = st.text_input("NCBI email, optional")
    ncbi_api_key = st.text_input("NCBI API key, optional", type="password")

    show_all = st.checkbox("Show all technical columns", value=False)

tab1, tab2, tab3 = st.tabs(["Run app", "Local database", "Template"])

with tab1:
    col1, col2 = st.columns(2)
    with col1:
        topic = st.text_input("Relation topic for literature search", placeholder="preeclampsia, gestational hypertension, preterm birth")
    with col2:
        input_mode = st.radio("Input mode", ["Single metabolite", "List of metabolites"], horizontal=True)

    if input_mode == "Single metabolite":
        met_text = st.text_input("Metabolite name", placeholder="Cysteine")
    else:
        met_text = st.text_area("Metabolite list", placeholder="Cysteine\nArginine\nCarnitine")

    if st.button("Run annotation and literature search", type="primary"):
        if not met_text.strip():
            st.error("Please enter at least one metabolite.")
            st.stop()

        mets = split_metabolites(met_text)
        ann_rows = []
        lit_rows = []
        progress = st.progress(0)
        status = st.empty()

        for i, m in enumerate(mets):
            status.write(f"Processing: {m}")
            ann_rows.append(annotate_one(m, local_db, use_pubchem, use_kegg, use_chebi, use_classyfire))

            if topic.strip():
                if use_pubmed:
                    lit_rows.extend(pubmed_search(m, topic, max_papers, ncbi_email, ncbi_api_key))
                    time.sleep(0.34)
                if use_epmc:
                    lit_rows.extend(europe_pmc_search(m, topic, max_papers))
                    time.sleep(0.34)

            progress.progress((i + 1) / len(mets))

        status.write("Completed.")

        ann_df = pd.DataFrame(ann_rows).fillna("Not available").replace("", "Not available")
        lit_df = pd.DataFrame(lit_rows).fillna("Not available").replace("", "Not available") if lit_rows else pd.DataFrame()

        priority_cols = [
            "input_metabolite",
            "best_compound_class",
            "best_pathway",
            "best_superpathway",
            "best_physiological_compartment",
            "compartment_source",
            "fallback_confidence",
            "hmdb_id",
            "kegg_id",
            "pubchem_cid",
            "pubchem_formula",
            "pubchem_molecular_weight",
            "pubchem_inchikey",
            "chebi_short_form",
            "classyfire_superclass",
            "classyfire_class",
            "pubchem_status",
            "kegg_status",
            "chebi_status",
            "classyfire_status"
        ]

        st.subheader("1. Integrated annotation results")
        st.dataframe(ann_df if show_all else ann_df[[c for c in priority_cols if c in ann_df.columns]], use_container_width=True)

        st.download_button(
            "Download full annotation CSV",
            data=ann_df.to_csv(index=False).encode("utf-8"),
            file_name="integrated_metabolite_annotation_v4_no_blank.csv",
            mime="text/csv"
        )

        st.subheader("2. Literature results")
        if topic.strip():
            if lit_df.empty:
                st.info("No literature results found.")
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
                    if str(row.get("error", "")).strip() not in ["", "Not available", "nan"]:
                        st.error(row.get("error"))
                        continue
                    st.markdown(f"**{row.get('Title', 'Not available')}**")
                    st.write(f"{row.get('source', 'Not available')} | {row.get('Journal', 'Not available')} | {row.get('Year', 'Not available')}")
                    if row.get("URL", "Not available") != "Not available":
                        st.markdown(f"[Open record]({row.get('URL')})")
                    if row.get("Abstract", "Not available") != "Not available":
                        abstract = str(row.get("Abstract", ""))
                        st.write(abstract[:1200] + ("..." if len(abstract) > 1200 else ""))
                    st.divider()
        else:
            st.info("No relation topic entered, so literature search was skipped.")

with tab2:
    st.subheader("Current local curated database")
    st.dataframe(local_db, use_container_width=True)

with tab3:
    template = pd.DataFrame(columns=[
        "metabolite", "synonyms", "compound_class", "pathway", "superpathway",
        "physiological_compartment", "hmdb_id", "kegg_id", "pubchem_cid",
        "inchikey", "notes"
    ])
    st.subheader("Annotation CSV template")
    st.dataframe(template, use_container_width=True)
    st.download_button(
        "Download template CSV",
        data=template.to_csv(index=False).encode("utf-8"),
        file_name="metabolite_annotation_template.csv",
        mime="text/csv"
    )
