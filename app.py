
import streamlit as st
import pandas as pd
import requests
import xml.etree.ElementTree as ET
import re
import time
from urllib.parse import quote_plus

st.set_page_config(
    page_title="Metabolite Annotation v5 - Genes & Proteins",
    page_icon="🧬",
    layout="wide"
)

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
KEGG_BASE = "https://rest.kegg.jp"
NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OLS_BASE = "https://www.ebi.ac.uk/ols/api"
UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb/search"

HEADERS = {"User-Agent": "MetaboliteAnnotationGeneProteinV5/1.0"}

# ============================================================
# General helpers
# ============================================================
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

def no_blank_dict(row):
    for k in list(row.keys()):
        if is_blank(row[k]):
            row[k] = "Not available"
    return row

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

# ============================================================
# Local curated database
# ============================================================
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
                "notes": "Sulfur amino acid involved in glutathione synthesis and redox balance.",
                "associated_genes": "CBS; CTH; GCLC; GCLM; GSS; SLC7A11; GGT1",
                "associated_proteins": "Cystathionine beta-synthase; Cystathionine gamma-lyase; Glutamate-cysteine ligase catalytic subunit; Glutamate-cysteine ligase modifier subunit; Glutathione synthetase; Cystine/glutamate transporter; Gamma-glutamyltransferase 1"
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
                "notes": "Substrate for nitric oxide synthase and arginase.",
                "associated_genes": "NOS3; NOS2; NOS1; ARG1; ARG2; ASS1; ASL; SLC7A1",
                "associated_proteins": "Endothelial nitric oxide synthase; Inducible nitric oxide synthase; Neuronal nitric oxide synthase; Arginase-1; Arginase-2; Argininosuccinate synthase; Argininosuccinate lyase; Cationic amino acid transporter 1"
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
                "notes": "Transports long-chain fatty acids into mitochondria.",
                "associated_genes": "CPT1A; CPT1B; CPT2; SLC22A5; CACT/SLC25A20; CRAT",
                "associated_proteins": "Carnitine palmitoyltransferase 1A; Carnitine palmitoyltransferase 1B; Carnitine palmitoyltransferase 2; Organic cation/carnitine transporter 2; Carnitine/acylcarnitine translocase; Carnitine O-acetyltransferase"
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
                "notes": "Precursor of kynurenine, serotonin, melatonin and microbial indoles.",
                "associated_genes": "IDO1; IDO2; TDO2; TPH1; TPH2; KYNU; KMO; AANAT",
                "associated_proteins": "Indoleamine 2,3-dioxygenase 1; Indoleamine 2,3-dioxygenase 2; Tryptophan 2,3-dioxygenase; Tryptophan hydroxylase 1; Tryptophan hydroxylase 2; Kynureninase; Kynurenine 3-monooxygenase; Serotonin N-acetyltransferase"
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
                "notes": "Important for membrane phospholipids and methylation biology.",
                "associated_genes": "CHKA; CHKB; PCYT1A; PEMT; BHMT; SLC44A1",
                "associated_proteins": "Choline kinase alpha; Choline kinase beta; Choline-phosphate cytidylyltransferase A; Phosphatidylethanolamine N-methyltransferase; Betaine-homocysteine S-methyltransferase; Choline transporter-like protein 1"
            }
        ])

    required = [
        "metabolite", "synonyms", "compound_class", "pathway", "superpathway",
        "physiological_compartment", "hmdb_id", "kegg_id", "pubchem_cid",
        "inchikey", "notes", "associated_genes", "associated_proteins"
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

# ============================================================
# Online annotation databases
# ============================================================
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

# ============================================================
# Fallback biological annotation
# ============================================================
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

def no_blank_annotation_row(row):
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

    return no_blank_dict(row)

# ============================================================
# Annotation integration
# ============================================================
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
        "local_associated_genes": local.get("associated_genes", ""),
        "local_associated_proteins": local.get("associated_proteins", ""),

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
    return no_blank_annotation_row(row)

# ============================================================
# Gene/protein association logic
# ============================================================
def split_semicolon(x):
    if is_blank(x):
        return []
    return [i.strip() for i in re.split(r"[;,|]+", str(x)) if i.strip() and i.strip() != "Not available"]

def fallback_genes_for_metabolite(name, annotation_row=None):
    text = name.lower()
    if annotation_row:
        text += " " + str(annotation_row.get("best_pathway", "")).lower()
        text += " " + str(annotation_row.get("best_compound_class", "")).lower()

    gene_map = []

    def add(gene, protein, function, pathway, source="Rule-based curated fallback", confidence="Moderate"):
        gene_map.append({
            "gene_name": gene,
            "protein_name": protein,
            "protein_details": function,
            "associated_pathway_or_process": pathway,
            "association_source": source,
            "association_confidence": confidence
        })

    if "cysteine" in text or "glutathione" in text or "sulfur" in text:
        add("CBS", "Cystathionine beta-synthase", "Catalyzes conversion reactions in transsulfuration and cysteine-related sulfur amino acid metabolism.", "Cysteine and methionine metabolism; transsulfuration")
        add("CTH", "Cystathionine gamma-lyase", "Generates cysteine and hydrogen sulfide-related sulfur metabolites.", "Cysteine and methionine metabolism")
        add("GCLC", "Glutamate-cysteine ligase catalytic subunit", "Rate-limiting enzyme for glutathione synthesis using cysteine as a substrate.", "Glutathione metabolism; redox regulation")
        add("GCLM", "Glutamate-cysteine ligase modifier subunit", "Regulatory subunit of glutamate-cysteine ligase involved in glutathione biosynthesis.", "Glutathione metabolism; redox regulation")
        add("GSS", "Glutathione synthetase", "Catalyzes final step in glutathione biosynthesis.", "Glutathione metabolism")
        add("SLC7A11", "Cystine/glutamate transporter", "Imports cystine for intracellular cysteine and glutathione production.", "Amino acid transport; redox homeostasis")
        add("GGT1", "Gamma-glutamyltransferase 1", "Participates in extracellular glutathione breakdown and cysteine recycling.", "Glutathione turnover")

    elif "arginine" in text or "nitric oxide" in text or "urea cycle" in text:
        add("NOS3", "Endothelial nitric oxide synthase", "Converts arginine to nitric oxide in vascular/endothelial biology.", "Nitric oxide biosynthesis")
        add("NOS2", "Inducible nitric oxide synthase", "Produces nitric oxide during inflammatory and immune responses.", "Nitric oxide biosynthesis; inflammation")
        add("ARG1", "Arginase-1", "Converts arginine to ornithine and urea, competing with nitric oxide synthases.", "Urea cycle; arginine metabolism")
        add("ARG2", "Arginase-2", "Mitochondrial arginase involved in arginine and ornithine metabolism.", "Arginine metabolism")
        add("ASS1", "Argininosuccinate synthase", "Catalyzes citrulline-to-argininosuccinate step in arginine biosynthesis/urea cycle.", "Urea cycle")
        add("ASL", "Argininosuccinate lyase", "Generates arginine and fumarate from argininosuccinate.", "Urea cycle")
        add("SLC7A1", "Cationic amino acid transporter 1", "Transports arginine and related cationic amino acids.", "Amino acid transport")

    elif "carnitine" in text or "fatty acid oxidation" in text:
        add("CPT1A", "Carnitine palmitoyltransferase 1A", "Controls entry of long-chain fatty acids into mitochondrial beta-oxidation.", "Carnitine shuttle; fatty acid oxidation")
        add("CPT1B", "Carnitine palmitoyltransferase 1B", "Muscle-enriched carnitine palmitoyltransferase involved in fatty acid oxidation.", "Carnitine shuttle")
        add("CPT2", "Carnitine palmitoyltransferase 2", "Converts acylcarnitines back to acyl-CoA inside mitochondria.", "Mitochondrial beta-oxidation")
        add("SLC22A5", "Organic cation/carnitine transporter 2", "High-affinity carnitine transporter important for systemic carnitine homeostasis.", "Carnitine transport")
        add("SLC25A20", "Carnitine/acylcarnitine translocase", "Transports acylcarnitines across the mitochondrial inner membrane.", "Mitochondrial carnitine shuttle")
        add("CRAT", "Carnitine O-acetyltransferase", "Catalyzes reversible transfer of acetyl groups between acetyl-CoA and carnitine.", "Acetylcarnitine metabolism")

    elif "tryptophan" in text or "kynurenine" in text or "serotonin" in text:
        add("IDO1", "Indoleamine 2,3-dioxygenase 1", "Initiates tryptophan catabolism through the kynurenine pathway during immune activation.", "Kynurenine pathway")
        add("IDO2", "Indoleamine 2,3-dioxygenase 2", "Tryptophan-catabolizing enzyme related to immune-metabolic regulation.", "Kynurenine pathway")
        add("TDO2", "Tryptophan 2,3-dioxygenase", "Liver-enriched enzyme initiating tryptophan degradation to kynurenine.", "Kynurenine pathway")
        add("TPH1", "Tryptophan hydroxylase 1", "Rate-limiting enzyme for peripheral serotonin synthesis.", "Serotonin biosynthesis")
        add("TPH2", "Tryptophan hydroxylase 2", "Neuronal tryptophan hydroxylase involved in serotonin production.", "Serotonin biosynthesis")
        add("KMO", "Kynurenine 3-monooxygenase", "Controls branch-point metabolism in the kynurenine pathway.", "Kynurenine pathway")

    elif "choline" in text or "phosphatidylcholine" in text:
        add("CHKA", "Choline kinase alpha", "Phosphorylates choline in phosphatidylcholine biosynthesis.", "Kennedy pathway; phosphatidylcholine metabolism")
        add("CHKB", "Choline kinase beta", "Choline kinase involved in membrane phospholipid biosynthesis.", "Phosphatidylcholine metabolism")
        add("PCYT1A", "Choline-phosphate cytidylyltransferase A", "Rate-limiting enzyme in phosphatidylcholine synthesis.", "Kennedy pathway")
        add("PEMT", "Phosphatidylethanolamine N-methyltransferase", "Converts phosphatidylethanolamine to phosphatidylcholine in liver.", "Phospholipid methylation")
        add("BHMT", "Betaine-homocysteine S-methyltransferase", "Links choline-derived betaine to one-carbon metabolism.", "One-carbon metabolism")
        add("SLC44A1", "Choline transporter-like protein 1", "Transports choline for phospholipid and methyl donor metabolism.", "Choline transport")

    elif "glucose" in text or "glycolysis" in text:
        add("HK1", "Hexokinase-1", "Phosphorylates glucose to glucose-6-phosphate.", "Glycolysis")
        add("GCK", "Glucokinase", "Liver/pancreatic glucose sensor enzyme.", "Glucose metabolism")
        add("SLC2A1", "Glucose transporter 1", "Facilitates cellular glucose uptake.", "Glucose transport")
        add("SLC2A4", "Glucose transporter 4", "Insulin-responsive glucose transporter.", "Glucose transport")
        add("G6PD", "Glucose-6-phosphate dehydrogenase", "Controls entry into pentose phosphate pathway and NADPH generation.", "Pentose phosphate pathway")

    elif "lactate" in text:
        add("LDHA", "L-lactate dehydrogenase A chain", "Converts pyruvate to lactate under glycolytic conditions.", "Lactate metabolism")
        add("LDHB", "L-lactate dehydrogenase B chain", "Catalyzes lactate-pyruvate interconversion.", "Lactate metabolism")
        add("SLC16A1", "Monocarboxylate transporter 1", "Transports lactate and other monocarboxylates.", "Lactate transport")
        add("SLC16A3", "Monocarboxylate transporter 4", "Exports lactate from glycolytic cells.", "Lactate transport")

    elif "bile acid" in text or "cholate" in text:
        add("CYP7A1", "Cholesterol 7-alpha-monooxygenase", "Rate-limiting enzyme in bile acid synthesis.", "Bile acid biosynthesis")
        add("CYP8B1", "Sterol 12-alpha-hydroxylase", "Controls cholic acid synthesis branch.", "Bile acid biosynthesis")
        add("SLC10A1", "Sodium/bile acid cotransporter", "Hepatic bile acid uptake transporter.", "Bile acid transport")
        add("ABCB11", "Bile salt export pump", "Exports bile acids from hepatocytes into bile.", "Bile acid transport")

    elif "steroid" in text or "cholesterol" in text or "progesterone" in text or "estradiol" in text:
        add("CYP11A1", "Cholesterol side-chain cleavage enzyme", "Converts cholesterol to pregnenolone.", "Steroidogenesis")
        add("HSD3B1", "3 beta-hydroxysteroid dehydrogenase type 1", "Catalyzes steroid hormone biosynthesis steps.", "Steroid hormone metabolism")
        add("CYP19A1", "Aromatase", "Converts androgens to estrogens.", "Estrogen biosynthesis")
        add("STAR", "Steroidogenic acute regulatory protein", "Moves cholesterol into mitochondria for steroidogenesis.", "Steroidogenesis")

    else:
        add("Not specifically identified", "Not specifically identified", "No confident metabolite-specific gene/protein association was found. Treat as likely exogenous, xenobiotic, microbial, or unclassified until manually curated.", "Unclassified/exogenous metabolism", "Fallback classification", "Low")

    return gene_map

@st.cache_data(show_spinner=False, ttl=86400)
def uniprot_by_ec(ec_number, max_results=5):
    if is_blank(ec_number):
        return []
    query = f'(ec:{ec_number}) AND (organism_id:9606)'
    params = {
        "query": query,
        "format": "json",
        "size": max_results,
        "fields": "accession,id,protein_name,gene_names,organism_name,cc_function,ec"
    }
    r, err = safe_get(UNIPROT_BASE, params=params, timeout=40)
    if r is None:
        return []
    try:
        data = r.json()
        results = data.get("results", [])
    except Exception:
        return []

    rows = []
    for item in results:
        accession = item.get("primaryAccession", "")
        protein_desc = item.get("proteinDescription", {})
        protein_name = ""
        try:
            protein_name = protein_desc.get("recommendedName", {}).get("fullName", {}).get("value", "")
        except Exception:
            protein_name = ""

        genes = []
        for g in item.get("genes", []) or []:
            if isinstance(g, dict):
                gn = g.get("geneName", {}).get("value", "")
                if gn:
                    genes.append(gn)

        comments = item.get("comments", []) or []
        functions = []
        for c in comments:
            if c.get("commentType") == "FUNCTION":
                for t in c.get("texts", []) or []:
                    val = t.get("value", "")
                    if val:
                        functions.append(val)

        rows.append({
            "gene_name": "; ".join(genes) or "Not available",
            "protein_name": protein_name or "Not available",
            "uniprot_id": accession or "Not available",
            "ec_number": ec_number,
            "protein_details": " ".join(functions)[:1200] if functions else "Function details not available from UniProt record",
            "associated_pathway_or_process": f"EC-linked enzyme activity from KEGG compound enzyme field: {ec_number}",
            "association_source": "KEGG enzyme field → UniProt human protein lookup",
            "association_confidence": "High if EC mapping is correct for the metabolite pathway"
        })
    return rows

@st.cache_data(show_spinner=False, ttl=86400)
def uniprot_by_gene(gene_name):
    if is_blank(gene_name) or gene_name == "Not specifically identified":
        return None
    query = f'(gene_exact:{gene_name}) AND (organism_id:9606)'
    params = {
        "query": query,
        "format": "json",
        "size": 1,
        "fields": "accession,id,protein_name,gene_names,organism_name,cc_function,ec"
    }
    r, err = safe_get(UNIPROT_BASE, params=params, timeout=30)
    if r is None:
        return None
    try:
        results = r.json().get("results", [])
        if not results:
            return None
        item = results[0]
    except Exception:
        return None

    accession = item.get("primaryAccession", "")
    protein_name = ""
    try:
        protein_name = item.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
    except Exception:
        protein_name = ""

    ec_numbers = []
    try:
        rec = item.get("proteinDescription", {}).get("recommendedName", {})
        for ec in rec.get("ecNumbers", []) or []:
            if ec.get("value"):
                ec_numbers.append(ec.get("value"))
    except Exception:
        pass

    functions = []
    for c in item.get("comments", []) or []:
        if c.get("commentType") == "FUNCTION":
            for t in c.get("texts", []) or []:
                if t.get("value"):
                    functions.append(t.get("value"))

    return {
        "uniprot_id": accession or "Not available",
        "protein_name_uniprot": protein_name or "Not available",
        "ec_number_uniprot": "; ".join(ec_numbers) if ec_numbers else "Not available",
        "protein_details_uniprot": " ".join(functions)[:1200] if functions else "Not available",
        "uniprot_url": f"https://www.uniprot.org/uniprotkb/{accession}/entry" if accession else "Not available"
    }

def build_gene_protein_rows(metabolite, ann_row, enrich_uniprot=True):
    rows = []

    # 1. KEGG EC enzyme → UniProt
    ec_numbers = split_semicolon(ann_row.get("kegg_enzyme", ""))
    for ec in ec_numbers[:12]:
        for uni in uniprot_by_ec(ec, max_results=4):
            row = {
                "metabolite": metabolite,
                "gene_name": uni.get("gene_name", ""),
                "protein_name": uni.get("protein_name", ""),
                "uniprot_id": uni.get("uniprot_id", ""),
                "ec_number": uni.get("ec_number", ""),
                "protein_details": uni.get("protein_details", ""),
                "associated_pathway_or_process": uni.get("associated_pathway_or_process", ""),
                "association_source": uni.get("association_source", ""),
                "association_confidence": uni.get("association_confidence", ""),
                "uniprot_url": f"https://www.uniprot.org/uniprotkb/{uni.get('uniprot_id')}/entry" if not is_blank(uni.get("uniprot_id")) else "Not available"
            }
            rows.append(no_blank_dict(row))

    # 2. Local associated gene/protein lists
    local_genes = split_semicolon(ann_row.get("local_associated_genes", ""))
    local_proteins = split_semicolon(ann_row.get("local_associated_proteins", ""))
    for idx, gene in enumerate(local_genes):
        protein = local_proteins[idx] if idx < len(local_proteins) else "Not available"
        row = {
            "metabolite": metabolite,
            "gene_name": gene,
            "protein_name": protein,
            "uniprot_id": "Not available",
            "ec_number": "Not available",
            "protein_details": "Locally curated metabolite-associated gene/protein. Details can be enriched from UniProt when available.",
            "associated_pathway_or_process": first_value(ann_row.get("best_pathway", ""), ann_row.get("best_superpathway", "")),
            "association_source": "Local curated annotation table",
            "association_confidence": "High if curated for the study; verify before publication",
            "uniprot_url": "Not available"
        }
        if enrich_uniprot:
            uni = uniprot_by_gene(gene)
            if uni:
                row["uniprot_id"] = first_value(uni.get("uniprot_id", ""), row["uniprot_id"])
                row["protein_name"] = first_value(uni.get("protein_name_uniprot", ""), row["protein_name"])
                row["ec_number"] = first_value(uni.get("ec_number_uniprot", ""), row["ec_number"])
                row["protein_details"] = first_value(uni.get("protein_details_uniprot", ""), row["protein_details"])
                row["uniprot_url"] = first_value(uni.get("uniprot_url", ""), row["uniprot_url"])
        rows.append(no_blank_dict(row))

    # 3. Rule-based fallback genes
    fallback = fallback_genes_for_metabolite(metabolite, ann_row)
    for fb in fallback:
        gene = fb.get("gene_name", "")
        row = {
            "metabolite": metabolite,
            "gene_name": gene,
            "protein_name": fb.get("protein_name", ""),
            "uniprot_id": "Not available",
            "ec_number": "Not available",
            "protein_details": fb.get("protein_details", ""),
            "associated_pathway_or_process": fb.get("associated_pathway_or_process", ""),
            "association_source": fb.get("association_source", ""),
            "association_confidence": fb.get("association_confidence", ""),
            "uniprot_url": "Not available"
        }
        if enrich_uniprot and gene != "Not specifically identified":
            uni = uniprot_by_gene(gene)
            if uni:
                row["uniprot_id"] = first_value(uni.get("uniprot_id", ""), row["uniprot_id"])
                row["protein_name"] = first_value(uni.get("protein_name_uniprot", ""), row["protein_name"])
                row["ec_number"] = first_value(uni.get("ec_number_uniprot", ""), row["ec_number"])
                row["protein_details"] = first_value(uni.get("protein_details_uniprot", ""), row["protein_details"])
                row["uniprot_url"] = first_value(uni.get("uniprot_url", ""), row["uniprot_url"])
        rows.append(no_blank_dict(row))

    # Deduplicate by metabolite + gene + protein + source priority
    if not rows:
        rows = [{
            "metabolite": metabolite,
            "gene_name": "Not specifically identified",
            "protein_name": "Not specifically identified",
            "uniprot_id": "Not available",
            "ec_number": "Not available",
            "protein_details": "No confident gene/protein association was found. Treat as exogenous, xenobiotic, microbial, or unclassified until manually curated.",
            "associated_pathway_or_process": first_value(ann_row.get("best_pathway", ""), "Unclassified metabolism"),
            "association_source": "Fallback classification",
            "association_confidence": "Low",
            "uniprot_url": "Not available"
        }]

    df = pd.DataFrame(rows).fillna("Not available").replace("", "Not available")
    if not df.empty:
        df["dedup_key"] = df["metabolite"].astype(str) + "|" + df["gene_name"].astype(str) + "|" + df["protein_name"].astype(str)
        df = df.drop_duplicates(subset=["dedup_key"]).drop(columns=["dedup_key"])
    return df.to_dict(orient="records")

# ============================================================
# Literature search
# ============================================================
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

    fetch_params = {"db": "pubmed", "id": ",".join(pmids), "retmode": "xml"}
    if email:
        fetch_params["email"] = email
    if api_key:
        fetch_params["api_key"] = api_key

    fr, ferr = safe_get(f"{NCBI_EUTILS}/efetch.fcgi", params=fetch_params)
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

# ============================================================
# UI
# ============================================================
st.title("🧬 Metabolite Annotation v5")
st.caption("No-blank output + gene/protein association tab")

with st.sidebar:
    st.header("Settings")
    uploaded_db = st.file_uploader("Optional curated annotation CSV", type=["csv"])
    local_db = load_local_db(uploaded_db)

    st.subheader("Annotation sources")
    use_pubchem = st.checkbox("Search PubChem", value=True)
    use_kegg = st.checkbox("Search KEGG", value=True)
    use_chebi = st.checkbox("Search ChEBI via EBI OLS", value=True)
    use_classyfire = st.checkbox("Search ClassyFire", value=True)

    st.subheader("Gene/protein settings")
    build_gene_tab = st.checkbox("Build gene/protein associations", value=True)
    enrich_uniprot = st.checkbox("Enrich genes/proteins from UniProt", value=True)

    st.subheader("Literature sources")
    use_pubmed = st.checkbox("Search PubMed", value=True)
    use_epmc = st.checkbox("Search Europe PMC", value=True)
    max_papers = st.slider("Maximum papers per source", 3, 30, 10)

    st.subheader("NCBI")
    ncbi_email = st.text_input("NCBI email, optional")
    ncbi_api_key = st.text_input("NCBI API key, optional", type="password")

    show_all = st.checkbox("Show all technical annotation columns", value=False)

tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "Run app",
    "Gene/protein associations",
    "Literature results",
    "Local database",
    "Template"
])

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

    if st.button("Run annotation, gene/protein and literature search", type="primary"):
        if not met_text.strip():
            st.error("Please enter at least one metabolite.")
            st.stop()

        mets = split_metabolites(met_text)
        ann_rows = []
        lit_rows = []
        gene_rows = []

        progress = st.progress(0)
        status = st.empty()

        for i, m in enumerate(mets):
            status.write(f"Processing: {m}")
            ann = annotate_one(m, local_db, use_pubchem, use_kegg, use_chebi, use_classyfire)
            ann_rows.append(ann)

            if build_gene_tab:
                gene_rows.extend(build_gene_protein_rows(m, ann, enrich_uniprot=enrich_uniprot))

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
        gene_df = pd.DataFrame(gene_rows).fillna("Not available").replace("", "Not available") if gene_rows else pd.DataFrame()
        lit_df = pd.DataFrame(lit_rows).fillna("Not available").replace("", "Not available") if lit_rows else pd.DataFrame()

        st.session_state["ann_df"] = ann_df
        st.session_state["gene_df"] = gene_df
        st.session_state["lit_df"] = lit_df

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
            file_name="integrated_metabolite_annotation_v5.csv",
            mime="text/csv"
        )

        if build_gene_tab:
            st.success("Gene/protein associations were created. Open the Gene/protein associations tab.")

        if topic.strip():
            st.success("Literature results were created. Open the Literature results tab.")
        else:
            st.info("No relation topic entered, so literature search was skipped.")

with tab2:
    st.subheader("Gene/protein associations")
    st.markdown(
        """
This tab links metabolites to genes/proteins using KEGG enzyme numbers, UniProt human protein records,
local curated annotations and rule-based fallback associations. Key fields are not left blank.
"""
    )

    gene_df = st.session_state.get("gene_df", pd.DataFrame())
    if gene_df.empty:
        st.info("Run the app first with 'Build gene/protein associations' enabled.")
    else:
        st.dataframe(gene_df, use_container_width=True)
        st.download_button(
            "Download gene/protein associations CSV",
            data=gene_df.to_csv(index=False).encode("utf-8"),
            file_name="metabolite_gene_protein_associations_v5.csv",
            mime="text/csv"
        )

        st.subheader("Gene/protein detail cards")
        for _, row in gene_df.head(50).iterrows():
            st.markdown(f"**{row.get('gene_name', 'Not available')} — {row.get('protein_name', 'Not available')}**")
            st.write(f"Metabolite: {row.get('metabolite', 'Not available')}")
            st.write(f"UniProt: {row.get('uniprot_id', 'Not available')} | EC: {row.get('ec_number', 'Not available')}")
            st.write(f"Process/pathway: {row.get('associated_pathway_or_process', 'Not available')}")
            st.write(f"Details: {row.get('protein_details', 'Not available')}")
            st.caption(f"Source: {row.get('association_source', 'Not available')} | Confidence: {row.get('association_confidence', 'Not available')}")
            if row.get("uniprot_url", "Not available") != "Not available":
                st.markdown(f"[Open UniProt record]({row.get('uniprot_url')})")
            st.divider()

with tab3:
    st.subheader("Literature results")
    lit_df = st.session_state.get("lit_df", pd.DataFrame())
    if lit_df.empty:
        st.info("No literature results available. Run the app with a relation topic.")
    else:
        st.dataframe(lit_df, use_container_width=True)
        st.download_button(
            "Download literature CSV",
            data=lit_df.to_csv(index=False).encode("utf-8"),
            file_name="metabolite_literature_results_v5.csv",
            mime="text/csv"
        )

        st.subheader("Paper cards")
        for _, row in lit_df.head(50).iterrows():
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

with tab4:
    st.subheader("Current local curated database")
    st.dataframe(local_db, use_container_width=True)

with tab5:
    template = pd.DataFrame(columns=[
        "metabolite", "synonyms", "compound_class", "pathway", "superpathway",
        "physiological_compartment", "hmdb_id", "kegg_id", "pubchem_cid",
        "inchikey", "notes", "associated_genes", "associated_proteins"
    ])
    st.subheader("Annotation CSV template")
    st.dataframe(template, use_container_width=True)
    st.download_button(
        "Download template CSV",
        data=template.to_csv(index=False).encode("utf-8"),
        file_name="metabolite_annotation_template_v5.csv",
        mime="text/csv"
    )
