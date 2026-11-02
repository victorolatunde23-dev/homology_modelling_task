"""
Homology Template Finder
========================
Takes a protein (UniProt accession, PDB ID[_chain], or raw/FASTA sequence),
finds all homologs with experimentally solved structures in the PDB, computes
template-selection metrics, ranks them with a transparent composite score and
writes everything to a formatted Excel workbook.

Designed for Google Colab (needs: requests, pandas, numpy, openpyxl;
Biopython is optional and only used for the NCBI BLAST fallback / realignment).
"""

import io
import json
import math
import os
import re
import time
import zipfile
from collections import defaultdict
from datetime import datetime

import numpy as np
import pandas as pd
import requests

# ---------------------------------------------------------------------------
# Default configuration (override any key when calling run_pipeline)
# ---------------------------------------------------------------------------
DEFAULT_CONFIG = {
    "query": "P69905",              # UniProt acc, PDB ID (e.g. 4HHB or 4HHB_A) or sequence/FASTA
    "search_engine": "rcsb",        # "rcsb" (fast MMseqs2) or "ncbi" (BLASTp vs pdb, slower)
    "evalue_cutoff": 10.0,          # loose -> finds remote homologs; tighten to 1e-3 for confident ones
    "min_identity_pct": 0,          # pre-filter on sequence identity (%)
    "max_identity_pct": 100,        # set e.g. 99 to drop the query's own structures
    "min_query_coverage_pct": 0,    # pre-filter on query coverage (%)
    "max_hits": 5000,               # safety cap
    "exclude_pdb_ids": [],          # e.g. ["4HHB"]
    "exclude_query_entry": False,   # if the query is a PDB ID, drop that entry from results
    "prefer_organism": "",          # optional: "Homo sapiens" -> adds a 'same organism' flag
    "weights": {                    # composite score weights (normalised automatically)
        "identity": 0.30,
        "similarity": 0.05,
        "coverage": 0.25,
        "resolution": 0.15,
        "completeness": 0.10,
        "r_free": 0.075,
        "geometry": 0.075,
    },
    "n_alignments_in_excel": 50,    # how many top hits get full alignments in the workbook
    "n_coverage_map": 50,           # how many top unique templates in the coverage map
    "coverage_bin": 10,             # residues per coverage-map column
    "download_top_n": 0,            # download top-N unique template structures (PDB/mmCIF) as a zip
    "output_prefix": "homology_templates",
}

RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_GRAPHQL = "https://data.rcsb.org/graphql"
HEADERS = {"User-Agent": "HomologyTemplateFinder/1.0 (Colab)"}

# Crystallisation additives / buffers / ions that are rarely relevant for modelling
COMMON_ADDITIVES = set("""
HOH DOD SO4 PO4 GOL EDO PEG PGE PG4 1PE P6G 2PE MPD ACT ACE FMT CL NA K BR IOD
MG CA ZN MN NI CO CU CD FE FE2 HG SR CS RB LI TRS EPE MES HEPES CIT TLA BME DTT
IMD NO3 NH4 SCN BU3 DMS EOH MOH IPA PEU UNX UNL BCN MLI SIN MAL TAR AZI 12P
15P P33 PE4 PE5 PE8 7PE XPE OXY CO3 BCT ACY FLC MRD
""".split())
# Metal ions are sometimes functionally important; list them separately
METAL_IONS = {"ZN", "MG", "CA", "MN", "FE", "FE2", "CU", "CO", "NI", "CD", "HG"}

# ---------------------------------------------------------------------------
# BLOSUM62 (embedded so Biopython is not required)
# ---------------------------------------------------------------------------
_B62 = """
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
"""


def _load_blosum62():
    lines = [l for l in _B62.strip("\n").split("\n") if l.strip()]
    cols = lines[0].split()
    m = {}
    for line in lines[1:]:
        parts = line.split()
        r = parts[0]
        for c, v in zip(cols, parts[1:]):
            m[(r, c)] = int(v)
    return m


BLOSUM62 = _load_blosum62()


def blosum(a, b):
    return BLOSUM62.get((a, b), BLOSUM62.get(("X", "X")))


# ---------------------------------------------------------------------------
# HTTP helpers
# ---------------------------------------------------------------------------
def _request(method, url, retries=4, timeout=120, **kw):
    kw.setdefault("headers", HEADERS)
    last = None
    for attempt in range(retries):
        try:
            r = requests.request(method, url, timeout=timeout, **kw)
            if r.status_code in (429, 500, 502, 503, 504):
                raise requests.HTTPError(f"HTTP {r.status_code}")
            return r
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"Request to {url} failed after {retries} attempts: {last}")


# ---------------------------------------------------------------------------
# 1. Query parsing
# ---------------------------------------------------------------------------
AA = set("ACDEFGHIKLMNPQRSTVWYBZXUO")
UNIPROT_RE = re.compile(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-\d+)?$")
PDB_RE = re.compile(r"^([0-9][A-Za-z0-9]{3})(?:[_:.\-]([A-Za-z0-9]{1,4}))?$")


def parse_fasta(text):
    name, seq = "query", []
    for line in text.strip().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if seq:
                break  # only the first record
            name = line[1:].strip() or "query"
        else:
            seq.append(re.sub(r"[^A-Za-z]", "", line))
    return name, "".join(seq).upper()


def fetch_query_sequence(query):
    """Return dict(name, sequence, source, pdb_entry, pdb_chain)."""
    q = query.strip()
    info = {"pdb_entry": None, "pdb_chain": None}

    if q.startswith(">") or "\n" in q:
        name, seq = parse_fasta(q)
        info.update(name=name, sequence=seq, source="User sequence")
        return info

    compact = re.sub(r"\s", "", q).upper()
    if len(compact) >= 20 and set(compact) <= AA and not UNIPROT_RE.match(compact):
        info.update(name="query", sequence=compact, source="User sequence")
        return info

    if UNIPROT_RE.match(compact):
        r = _request("GET", f"https://rest.uniprot.org/uniprotkb/{compact}.fasta")
        if r.status_code != 200 or not r.text.startswith(">"):
            raise ValueError(f"UniProt accession {compact} not found (HTTP {r.status_code}).")
        name, seq = parse_fasta(r.text)
        info.update(name=name, sequence=seq, source=f"UniProt {compact}")
        return info

    m = PDB_RE.match(q)
    if m:
        pdb_id, chain = m.group(1).upper(), m.group(2)
        r = _request("GET", f"https://www.rcsb.org/fasta/entry/{pdb_id}")
        if r.status_code != 200 or not r.text.startswith(">"):
            raise ValueError(f"PDB entry {pdb_id} not found (HTTP {r.status_code}).")
        records = []
        for block in r.text.strip().split(">")[1:]:
            lines = block.strip().splitlines()
            header, seq = lines[0], "".join(lines[1:]).upper()
            parts = header.split("|")
            chains_field = parts[1] if len(parts) > 1 else ""
            chains_field = re.sub(r"^Chains?\s+", "", chains_field)
            chains = []
            for c in chains_field.split(","):
                c = c.strip()
                mm = re.match(r"(\S+?)\[auth (\S+)\]", c)
                if mm:
                    chains += [mm.group(1), mm.group(2)]
                elif c:
                    chains.append(c)
            records.append({"header": header, "seq": seq, "chains": chains})
        prot = [rec for rec in records if set(rec["seq"]) <= AA and len(rec["seq"]) >= 10]
        if not prot:
            raise ValueError(f"No protein chains found in {pdb_id}.")
        chosen = prot[0]
        if chain:
            hit = [rec for rec in prot if chain in rec["chains"]]
            if not hit:
                raise ValueError(f"Chain {chain} not found in {pdb_id}. Chains: "
                                 + ", ".join(sorted({c for rec in prot for c in rec['chains']})))
            chosen = hit[0]
        info.update(name=chosen["header"], sequence=chosen["seq"], source=f"PDB {pdb_id}",
                    pdb_entry=pdb_id, pdb_chain=chain or (chosen["chains"][0] if chosen["chains"] else None))
        return info

    raise ValueError("Could not interpret the query. Give a UniProt accession (P69905), "
                     "a PDB ID (4HHB or 4HHB_A) or an amino-acid sequence / FASTA.")


# ---------------------------------------------------------------------------
# 2. Homolog search
# ---------------------------------------------------------------------------
def search_rcsb(sequence, evalue_cutoff=10.0, identity_cutoff=0.0, max_hits=5000):
    """MMseqs2 search of the PDB via the RCSB Search API. Returns list of hit dicts."""
    payload = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": evalue_cutoff,
                "identity_cutoff": identity_cutoff,
                "sequence_type": "protein",
                "value": sequence,
            },
        },
        "request_options": {
            "scoring_strategy": "sequence",
            "results_verbosity": "verbose",
            "paginate": {"start": 0, "rows": max_hits},
        },
        "return_type": "polymer_entity",
    }
    r = _request("POST", RCSB_SEARCH, json=payload)
    if r.status_code == 204:
        return []
    if r.status_code != 200:
        raise RuntimeError(f"RCSB search failed: HTTP {r.status_code}: {r.text[:500]}")
    data = r.json()
    hits = []
    for res in data.get("result_set", []):
        ident = res["identifier"]  # e.g. 4HHB_1
        for svc in res.get("services", []):
            for node in svc.get("nodes", []):
                for mc in node.get("match_context", []):
                    hits.append({
                        "entity_id": ident,
                        "pdb_id": ident.split("_")[0].upper(),
                        "evalue": mc.get("evalue"),
                        "bitscore": mc.get("bitscore"),
                        "search_identity": mc.get("sequence_identity"),
                        "query_beg": mc.get("query_beg"),
                        "query_end": mc.get("query_end"),
                        "subject_beg": mc.get("subject_beg"),
                        "subject_end": mc.get("subject_end"),
                        "subject_length": mc.get("subject_length"),
                        "query_aligned_seq": mc.get("query_aligned_seq"),
                        "subject_aligned_seq": mc.get("subject_aligned_seq"),
                        "hit_chains": None,
                    })
    # keep best HSP per entity
    best = {}
    for h in hits:
        k = h["entity_id"]
        if k not in best or (h["bitscore"] or 0) > (best[k]["bitscore"] or 0):
            best[k] = h
    return list(best.values())


def search_ncbi(sequence, evalue_cutoff=10.0, max_hits=5000):
    """BLASTp against NCBI 'pdb' database (slower, needs Biopython)."""
    from Bio.Blast import NCBIWWW, NCBIXML  # noqa: WPS433
    print("  Submitting BLASTp to NCBI (this can take a few minutes)...")
    handle = NCBIWWW.qblast("blastp", "pdb", sequence, expect=evalue_cutoff,
                            hitlist_size=min(max_hits, 5000), alignments=min(max_hits, 5000),
                            descriptions=min(max_hits, 5000))
    rec = NCBIXML.read(handle)
    raw = []
    for aln in rec.alignments:
        chains = re.findall(r"pdb\|(\w{4})\|(\w+)", aln.hit_id + " " + aln.hit_def)
        if not chains:
            continue
        hsp = max(aln.hsps, key=lambda h: h.bits)
        for pdb_id, chain in chains:
            raw.append({
                "pdb_id": pdb_id.upper(), "chain": chain,
                "evalue": hsp.expect, "bitscore": hsp.bits, "search_identity": None,
                "query_beg": hsp.query_start, "query_end": hsp.query_end,
                "subject_beg": hsp.sbjct_start, "subject_end": hsp.sbjct_end,
                "subject_length": aln.length,
                "query_aligned_seq": hsp.query, "subject_aligned_seq": hsp.sbjct,
            })
    if not raw:
        return []
    # map auth chains -> polymer entity ids
    pdb_ids = sorted({h["pdb_id"] for h in raw})
    chain2ent = {}
    q = """query($ids:[String!]!){ entries(entry_ids:$ids){ rcsb_id
           polymer_entities{ rcsb_id rcsb_polymer_entity_container_identifiers{ auth_asym_ids asym_ids } } } }"""
    for i in range(0, len(pdb_ids), 200):
        d = _graphql(q, {"ids": pdb_ids[i:i + 200]})
        for e in (d or {}).get("entries") or []:
            for pe in e.get("polymer_entities") or []:
                ci = pe.get("rcsb_polymer_entity_container_identifiers") or {}
                for c in (ci.get("auth_asym_ids") or []) + (ci.get("asym_ids") or []):
                    chain2ent.setdefault((e["rcsb_id"].upper(), c), pe["rcsb_id"])
    best = {}
    for h in raw:
        ent = chain2ent.get((h["pdb_id"], h["chain"]))
        if not ent:
            continue
        h = dict(h, entity_id=ent, hit_chains=h["chain"])
        h.pop("chain")
        if ent not in best or h["bitscore"] > best[ent]["bitscore"]:
            best[ent] = h
    return list(best.values())


# ---------------------------------------------------------------------------
# 3. Metadata from RCSB GraphQL
# ---------------------------------------------------------------------------
def _graphql(query, variables):
    r = _request("POST", RCSB_GRAPHQL, json={"query": query, "variables": variables})
    if r.status_code != 200:
        raise RuntimeError(f"GraphQL HTTP {r.status_code}: {r.text[:300]}")
    js = r.json()
    if js.get("errors") and not js.get("data"):
        raise RuntimeError("GraphQL error: " + json.dumps(js["errors"])[:500])
    if js.get("errors"):
        # partial data with field errors -> treat as schema problem
        raise RuntimeError("GraphQL field error: " + json.dumps(js["errors"])[:500])
    return js.get("data")


_ENTITY_CORE = """
    rcsb_id
    rcsb_polymer_entity { pdbx_description }
    entity_poly { pdbx_seq_one_letter_code_can rcsb_sample_sequence_length rcsb_mutation_count }
    rcsb_entity_source_organism { scientific_name ncbi_taxonomy_id }
    rcsb_polymer_entity_container_identifiers { auth_asym_ids uniprot_ids }
"""
_INSTANCES = """
    polymer_entity_instances {
      rcsb_polymer_entity_instance_container_identifiers { auth_asym_id asym_id }
      rcsb_polymer_instance_feature { type feature_positions { beg_seq_id end_seq_id } }
    }
"""
_ENTRY_CORE = """
      rcsb_id
      struct { title }
      exptl { method }
      rcsb_entry_info { resolution_combined polymer_entity_count_protein
                        deposited_polymer_entity_instance_count nonpolymer_bound_components }
      refine { ls_R_factor_R_free ls_R_factor_R_work }
      rcsb_accession_info { deposit_date initial_release_date }
      nonpolymer_entities { nonpolymer_comp { chem_comp { id name } } }
"""
_ENTRY_VALIDATION = """
      pdbx_vrpt_summary_geometry { clashscore percent_ramachandran_outliers percent_rotamer_outliers }
"""
_ENTRY_VALIDATION_OLD = """
      pdbx_vrpt_summary { clashscore percent_ramachandran_outliers percent_rotamer_outliers }
"""
_ENTRY_ASSEMBLY = """
      assemblies { rcsb_struct_symmetry { oligomeric_state stoichiometry kind } }
"""


def _build_entity_query(instances=True, validation="new", assembly=True):
    entry = _ENTRY_CORE
    if validation == "new":
        entry += _ENTRY_VALIDATION
    elif validation == "old":
        entry += _ENTRY_VALIDATION_OLD
    if assembly:
        entry += _ENTRY_ASSEMBLY
    return ("query($ids:[String!]!){ polymer_entities(entity_ids:$ids){"
            + _ENTITY_CORE + (_INSTANCES if instances else "")
            + " entry {" + entry + "} } }")


_QUERY_VARIANTS = [
    dict(instances=True, validation="new", assembly=True),
    dict(instances=True, validation="old", assembly=True),
    dict(instances=True, validation=None, assembly=True),
    dict(instances=False, validation=None, assembly=False),
]


def fetch_entity_metadata(entity_ids, batch=100, verbose=True):
    """Return {entity_id: metadata dict}."""
    out = {}
    variant_idx = 0
    ids = list(entity_ids)
    i = 0
    while i < len(ids):
        chunk = ids[i:i + batch]
        q = _build_entity_query(**_QUERY_VARIANTS[variant_idx])
        try:
            data = _graphql(q, {"ids": chunk})
        except RuntimeError as e:
            if "field" in str(e).lower() or "error" in str(e).lower():
                if variant_idx < len(_QUERY_VARIANTS) - 1:
                    variant_idx += 1
                    if verbose:
                        print(f"  GraphQL schema mismatch, retrying with reduced query (variant {variant_idx}).")
                    continue
            raise
        for pe in (data or {}).get("polymer_entities") or []:
            out[pe["rcsb_id"]] = pe
        i += batch
        if verbose:
            print(f"  metadata: {min(i, len(ids))}/{len(ids)} entities")
    return out


def _first(x, default=None):
    if isinstance(x, list):
        return x[0] if x else default
    return x if x is not None else default


def flatten_metadata(pe):
    """Turn a GraphQL polymer_entity record into flat fields."""
    entry = pe.get("entry") or {}
    info = entry.get("rcsb_entry_info") or {}
    refine = _first(entry.get("refine"), {}) or {}
    exptl = entry.get("exptl") or []
    methods = sorted({e.get("method") for e in exptl if e.get("method")})
    val = _first(entry.get("pdbx_vrpt_summary_geometry"), None) or entry.get("pdbx_vrpt_summary") or {}
    if isinstance(val, list):
        val = _first(val, {}) or {}
    org = _first(pe.get("rcsb_entity_source_organism"), {}) or {}
    ci = pe.get("rcsb_polymer_entity_container_identifiers") or {}
    poly = pe.get("entity_poly") or {}
    acc = entry.get("rcsb_accession_info") or {}

    ligands, metals = [], []
    for npe in entry.get("nonpolymer_entities") or []:
        cc = ((npe or {}).get("nonpolymer_comp") or {}).get("chem_comp") or {}
        cid = cc.get("id")
        if not cid:
            continue
        if cid in METAL_IONS:
            metals.append(cid)
        elif cid not in COMMON_ADDITIVES:
            ligands.append(f"{cid} ({cc.get('name', '')[:40]})")

    olig, stoich = None, None
    for asm in entry.get("assemblies") or []:
        for sym in (asm or {}).get("rcsb_struct_symmetry") or []:
            if sym and sym.get("kind") in (None, "Global Symmetry"):
                olig = olig or sym.get("oligomeric_state")
                st = sym.get("stoichiometry")
                stoich = stoich or (",".join(st) if isinstance(st, list) else st)
        if olig:
            break

    # unobserved residues per chain (entity numbering)
    chain_missing = {}
    for inst in pe.get("polymer_entity_instances") or []:
        ici = inst.get("rcsb_polymer_entity_instance_container_identifiers") or {}
        ch = ici.get("auth_asym_id") or ici.get("asym_id")
        segs = []
        for f in inst.get("rcsb_polymer_instance_feature") or []:
            if f.get("type") == "UNOBSERVED_RESIDUE_XYZ":
                for p in f.get("feature_positions") or []:
                    if p.get("beg_seq_id") is not None:
                        segs.append((p["beg_seq_id"], p.get("end_seq_id") or p["beg_seq_id"]))
        chain_missing[ch] = segs

    res = _first(info.get("resolution_combined"))
    return {
        "title": (entry.get("struct") or {}).get("title"),
        "description": (pe.get("rcsb_polymer_entity") or {}).get("pdbx_description"),
        "organism": org.get("scientific_name"),
        "taxid": org.get("ncbi_taxonomy_id"),
        "uniprot_ids": ";".join(ci.get("uniprot_ids") or []),
        "chains": ",".join(ci.get("auth_asym_ids") or []),
        "entity_sequence": (poly.get("pdbx_seq_one_letter_code_can") or "").replace("\n", ""),
        "template_length": poly.get("rcsb_sample_sequence_length"),
        "mutations": poly.get("rcsb_mutation_count"),
        "method": "; ".join(methods),
        "resolution": res,
        "r_free": refine.get("ls_R_factor_R_free"),
        "r_work": refine.get("ls_R_factor_R_work"),
        "clashscore": val.get("clashscore"),
        "rama_outliers_pct": val.get("percent_ramachandran_outliers"),
        "rotamer_outliers_pct": val.get("percent_rotamer_outliers"),
        "deposit_date": (acc.get("deposit_date") or "")[:10],
        "release_date": (acc.get("initial_release_date") or "")[:10],
        "n_protein_entities": info.get("polymer_entity_count_protein"),
        "n_chains_entry": info.get("deposited_polymer_entity_instance_count"),
        "oligomeric_state": olig,
        "stoichiometry": stoich,
        "ligands": "; ".join(sorted(set(ligands))),
        "metal_ions": ";".join(sorted(set(metals))),
        "_chain_missing": chain_missing,
    }


# ---------------------------------------------------------------------------
# 4. Alignment metrics
# ---------------------------------------------------------------------------
def local_align(a, b, gap_open=-11, gap_ext=-1):
    """Smith-Waterman with affine gaps (BLOSUM62). Uses Biopython if available."""
    try:
        from Bio import Align
        from Bio.Align import substitution_matrices
        al = Align.PairwiseAligner()
        al.mode = "local"
        al.substitution_matrix = substitution_matrices.load("BLOSUM62")
        al.open_gap_score, al.extend_gap_score = gap_open, gap_ext
        aln = al.align(a, b)[0]
        qa, sa = [], []
        (qs, _), (ss, _) = aln.aligned[0][0], aln.aligned[1][0]
        qi, si = qs, ss
        for (q0, q1), (s0, s1) in zip(*aln.aligned):
            if q0 > qi:
                qa.append(a[qi:q0]); sa.append("-" * (q0 - qi))
            if s0 > si:
                qa.append("-" * (s0 - si)); sa.append(b[si:s0])
            qa.append(a[q0:q1]); sa.append(b[s0:s1])
            qi, si = q1, s1
        return "".join(qa), "".join(sa), qs + 1, qi, ss + 1, si
    except ImportError:
        pass
    n, m = len(a), len(b)
    NEG = -10 ** 9
    H = np.zeros((n + 1, m + 1), dtype=np.int32)
    E = np.full((n + 1, m + 1), NEG, dtype=np.int64)
    F = np.full((n + 1, m + 1), NEG, dtype=np.int64)
    T = np.zeros((n + 1, m + 1), dtype=np.int8)
    best, bi, bj = 0, 0, 0
    for i in range(1, n + 1):
        ai = a[i - 1]
        for j in range(1, m + 1):
            E[i, j] = max(E[i, j - 1] + gap_ext, H[i, j - 1] + gap_open)
            F[i, j] = max(F[i - 1, j] + gap_ext, H[i - 1, j] + gap_open)
            d = H[i - 1, j - 1] + blosum(ai, b[j - 1])
            v = max(0, d, E[i, j], F[i, j])
            H[i, j] = v
            T[i, j] = 0 if v == 0 else (1 if v == d else (2 if v == F[i, j] else 3))
            if v > best:
                best, bi, bj = v, i, j
    i, j, qa, sa = bi, bj, [], []
    while i > 0 and j > 0 and H[i, j] > 0:
        t = T[i, j]
        if t == 1:
            qa.append(a[i - 1]); sa.append(b[j - 1]); i -= 1; j -= 1
        elif t == 2:
            qa.append(a[i - 1]); sa.append("-"); i -= 1
        else:
            qa.append("-"); sa.append(b[j - 1]); j -= 1
    return "".join(reversed(qa)), "".join(reversed(sa)), i + 1, bi, j + 1, bj


def alignment_stats(qa, sa):
    ident = pos = gaps_q = gaps_s = opens = 0
    prev_gap = False
    for x, y in zip(qa, sa):
        if x == "-" or y == "-":
            if x == "-":
                gaps_q += 1
            else:
                gaps_s += 1
            if not prev_gap:
                opens += 1
            prev_gap = True
            continue
        prev_gap = False
        if x == y:
            ident += 1
            pos += 1
        elif blosum(x, y) > 0:
            pos += 1
    L = len(qa)
    return {"aln_len": L, "identities": ident, "positives": pos,
            "gaps_in_query": gaps_q, "gaps_in_template": gaps_s, "gap_opens": opens}


def match_line(qa, sa):
    return "".join("|" if x == y and x != "-" else (":" if x != "-" and y != "-" and blosum(x, y) > 0
                                                    else " ") for x, y in zip(qa, sa))


def hssp_threshold(L):
    """Rost (1999) twilight-zone curve: identity (%) needed for a reliable homology at length L."""
    if L <= 0:
        return 100.0
    L = min(L, 450)
    return min(100.0, 480.0 * L ** (-0.32 * (1 + math.exp(-L / 1000.0))))


def query_positions_covered(qa, sa, qbeg):
    """Set of query residue indices (1-based) aligned to a template residue."""
    pos, qi = set(), qbeg
    for x, y in zip(qa, sa):
        if x != "-":
            if y != "-":
                pos.add(qi)
            qi += 1
    return pos


def template_positions_used(qa, sa, sbeg):
    pos, si = [], sbeg
    for x, y in zip(qa, sa):
        if y != "-":
            if x != "-":
                pos.append(si)
            si += 1
    return pos


# ---------------------------------------------------------------------------
# 5. Scoring
# ---------------------------------------------------------------------------
def _lin(x, lo, hi, invert=False):
    """Linear 0-1 scaling clipped; invert=True means lower is better."""
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return None
    v = (x - lo) / (hi - lo)
    v = min(1.0, max(0.0, v))
    return 1 - v if invert else v


def component_scores(row):
    meth = (row.get("method") or "").upper()
    s = {}
    s["identity"] = _lin(row["identity_pct"], 20, 100)
    s["similarity"] = _lin(row["similarity_pct"], 30, 100)
    s["coverage"] = _lin(row["query_coverage_pct"], 0, 100)
    res = row.get("resolution")
    if "NMR" in meth and (res is None or pd.isna(res)):
        s["resolution"] = 0.35
    else:
        r = _lin(row.get("resolution"), 1.5, 3.5, invert=True)
        s["resolution"] = 0.3 if r is None else r
    c = row.get("aligned_region_observed_pct")
    s["completeness"] = 0.5 if c is None or pd.isna(c) else _lin(c, 70, 100)
    rf = _lin(row.get("r_free"), 0.20, 0.32, invert=True)
    s["r_free"] = 0.5 if rf is None else rf
    g = [v for v in (_lin(row.get("clashscore"), 5, 40, invert=True),
                     _lin(row.get("rama_outliers_pct"), 0.5, 5, invert=True)) if v is not None]
    s["geometry"] = float(np.mean(g)) if g else 0.5
    return s


def confidence_tier(identity, coverage):
    if identity >= 50 and coverage >= 70:
        return "A - High (expect ~1-2 Å Cα RMSD core)"
    if identity >= 30 and coverage >= 50:
        return "B - Good (standard homology modelling)"
    if identity >= 20:
        return "C - Twilight zone (check alignment / use multiple templates)"
    return "D - Remote (fold recognition / threading only)"


def safe(v):
    return None if v is None or (isinstance(v, float) and math.isnan(v)) else v


# ---------------------------------------------------------------------------
# 6. Main pipeline
# ---------------------------------------------------------------------------
def build_table(query_info, hits, meta, cfg):
    qseq = query_info["sequence"]
    qlen = len(qseq)
    rows, coverage_sets, aln_store = [], {}, {}
    pref_org = (cfg.get("prefer_organism") or "").strip().lower()

    for h in hits:
        pe = meta.get(h["entity_id"])
        if pe is None:
            continue
        m = flatten_metadata(pe)
        qa, sa = h.get("query_aligned_seq"), h.get("subject_aligned_seq")
        if not qa or not sa or len(qa) != len(sa):
            if not m["entity_sequence"]:
                continue
            qa, sa, qb, qe, sb, se = local_align(qseq, m["entity_sequence"])
            h.update(query_beg=qb, query_end=qe, subject_beg=sb, subject_end=se)
        qa, sa = qa.upper(), sa.upper()
        st = alignment_stats(qa, sa)
        aligned_cols = st["aln_len"] - st["gaps_in_query"] - st["gaps_in_template"]
        qcov_set = query_positions_covered(qa, sa, h["query_beg"])
        tpos = template_positions_used(qa, sa, h["subject_beg"])

        # observed fraction of template residues used by the alignment (best chain)
        best_chain, best_obs = None, None
        if m["_chain_missing"] and tpos:
            for ch, segs in m["_chain_missing"].items():
                miss = sum(1 for p in tpos if any(b <= p <= e for b, e in segs))
                obs = 100.0 * (1 - miss / len(tpos))
                if best_obs is None or obs > best_obs:
                    best_chain, best_obs = ch, obs
        if h.get("hit_chains"):
            best_chain = best_chain or h["hit_chains"]

        identity_pct = 100.0 * st["identities"] / st["aln_len"] if st["aln_len"] else 0
        L_ungapped = aligned_cols
        hssp = hssp_threshold(L_ungapped)
        tlen = m["template_length"] or h.get("subject_length")
        row = {
            "entity_id": h["entity_id"],
            "pdb_id": h["pdb_id"],
            "best_chain": best_chain or (m["chains"].split(",")[0] if m["chains"] else None),
            "all_chains": m["chains"],
            "description": m["description"],
            "organism": m["organism"],
            "uniprot_ids": m["uniprot_ids"],
            "identity_pct": round(identity_pct, 2),
            "similarity_pct": round(100.0 * st["positives"] / st["aln_len"], 2) if st["aln_len"] else 0,
            "identity_vs_query_len_pct": round(100.0 * st["identities"] / qlen, 2),
            "query_coverage_pct": round(100.0 * len(qcov_set) / qlen, 2),
            "template_coverage_pct": round(100.0 * len(tpos) / tlen, 2) if tlen else None,
            "query_start": h["query_beg"],
            "query_end": h["query_end"],
            "template_start": h["subject_beg"],
            "template_end": h["subject_end"],
            "query_length": qlen,
            "template_length": tlen,
            "aligned_residues": L_ungapped,
            "alignment_length": st["aln_len"],
            "gap_opens": st["gap_opens"],
            "gaps_in_query": st["gaps_in_query"],
            "gaps_in_template": st["gaps_in_template"],
            "evalue": safe(h.get("evalue")),
            "bitscore": safe(h.get("bitscore")),
            "hssp_threshold_pct": round(hssp, 1),
            "hssp_distance": round(identity_pct - hssp, 1),
            "method": m["method"],
            "resolution": m["resolution"],
            "r_free": m["r_free"],
            "r_work": m["r_work"],
            "clashscore": m["clashscore"],
            "rama_outliers_pct": m["rama_outliers_pct"],
            "rotamer_outliers_pct": m["rotamer_outliers_pct"],
            "aligned_region_observed_pct": round(best_obs, 1) if best_obs is not None else None,
            "mutations": m["mutations"],
            "oligomeric_state": m["oligomeric_state"],
            "stoichiometry": m["stoichiometry"],
            "n_protein_entities": m["n_protein_entities"],
            "n_chains_entry": m["n_chains_entry"],
            "ligands": m["ligands"],
            "metal_ions": m["metal_ions"],
            "release_date": m["release_date"],
            "title": m["title"],
            "same_organism": (pref_org in (m["organism"] or "").lower()) if pref_org else None,
            "_entity_sequence": m["entity_sequence"],
        }
        rows.append(row)
        coverage_sets[h["entity_id"]] = qcov_set
        aln_store[h["entity_id"]] = (qa, sa)

    df = pd.DataFrame(rows)
    if df.empty:
        return df, coverage_sets, aln_store

    # filters
    df = df[(df.identity_pct >= cfg["min_identity_pct"]) & (df.identity_pct <= cfg["max_identity_pct"])
            & (df.query_coverage_pct >= cfg["min_query_coverage_pct"])]
    excl = {x.upper() for x in cfg.get("exclude_pdb_ids") or []}
    if cfg.get("exclude_query_entry") and query_info.get("pdb_entry"):
        excl.add(query_info["pdb_entry"].upper())
    if excl:
        df = df[~df.pdb_id.isin(excl)]

    # scoring
    w = cfg["weights"]
    wsum = sum(w.values())
    comps = df.apply(lambda r: component_scores(r.to_dict()), axis=1)
    for k in w:
        df[f"score_{k}"] = comps.apply(lambda s: round(s[k], 3))
    df["template_score"] = sum(df[f"score_{k}"] * (w[k] / wsum) for k in w) * 100
    df["template_score"] = df["template_score"].round(1)
    df["confidence_tier"] = [confidence_tier(i, c) for i, c in zip(df.identity_pct, df.query_coverage_pct)]
    df["homology_reliable"] = df["hssp_distance"] > 0

    warn = []
    for _, r in df.iterrows():
        f = []
        if r.query_coverage_pct < 50:
            f.append("low coverage")
        if r.identity_pct < 25:
            f.append("low identity")
        if r.resolution and r.resolution > 3.0:
            f.append("low resolution")
        if "NMR" in (r.method or "").upper():
            f.append("NMR model")
        if r.aligned_region_observed_pct is not None and r.aligned_region_observed_pct < 90:
            f.append("missing residues in aligned region")
        if r.mutations and r.mutations > 0:
            f.append(f"{int(r.mutations)} mutation(s)")
        if r.n_protein_entities and r.n_protein_entities > 1:
            f.append("complex with other proteins")
        if not r.homology_reliable:
            f.append("below HSSP curve")
        warn.append("; ".join(f))
    df["warnings"] = warn

    df = df.sort_values(["template_score", "bitscore"], ascending=[False, False]).reset_index(drop=True)
    df.insert(0, "rank", range(1, len(df) + 1))
    return df, coverage_sets, aln_store


def unique_templates(df):
    """One representative per identical template sequence (best score)."""
    if df.empty:
        return df
    g = df.groupby("_entity_sequence", sort=False)
    counts = g["pdb_id"].transform("count")
    allids = g["pdb_id"].transform(lambda s: ",".join(dict.fromkeys(s)))
    d = df.assign(n_structures_same_sequence=counts, other_pdb_ids_same_sequence=allids)
    u = d.drop_duplicates("_entity_sequence", keep="first").copy()
    u["other_pdb_ids_same_sequence"] = [
        ",".join([x for x in ids.split(",") if x != p][:40]) for ids, p in zip(u.other_pdb_ids_same_sequence, u.pdb_id)]
    u = u.reset_index(drop=True)
    u.insert(0, "unique_rank", range(1, len(u) + 1))
    return u


def uncovered_regions(sets, qlen):
    allpos = set().union(*sets) if sets else set()
    regs, start = [], None
    for i in range(1, qlen + 2):
        if i <= qlen and i not in allpos:
            start = start or i
        elif start:
            regs.append((start, i - 1))
            start = None
    return regs


def greedy_multitemplate(uniq, coverage_sets, qlen, max_templates=5, min_gain=10):
    """Pick a small set of complementary templates that maximise coverage (score-ordered greedy)."""
    chosen, covered = [], set()
    cands = uniq.head(200)
    for _ in range(max_templates):
        best, best_gain = None, 0
        for _, r in cands.iterrows():
            if r.entity_id in [c["entity_id"] for c in chosen]:
                continue
            gain = len(coverage_sets[r.entity_id] - covered)
            # favour quality: first template = top ranked
            if not chosen:
                best, best_gain = r, gain
                break
            if gain > best_gain:
                best, best_gain = r, gain
        if best is None or (chosen and best_gain < min_gain):
            break
        covered |= coverage_sets[best.entity_id]
        chosen.append({"entity_id": best.entity_id, "pdb_id": best.pdb_id, "chain": best.best_chain,
                       "query_range": f"{best.query_start}-{best.query_end}",
                       "identity_pct": best.identity_pct, "new_residues_covered": best_gain,
                       "cumulative_coverage_pct": round(100 * len(covered) / qlen, 1)})
    return pd.DataFrame(chosen)


# ---------------------------------------------------------------------------
# 7. Excel writer
# ---------------------------------------------------------------------------
COLUMN_DOCS = [
    ("rank", "Overall rank by template_score (ties broken by bit score)."),
    ("template_score", "Composite 0-100 score = weighted sum of the score_* components (weights on 'Settings' sheet)."),
    ("confidence_tier", "A: ≥50% id & ≥70% cov; B: ≥30% id & ≥50% cov; C: 20-30% id (twilight zone); D: <20% id."),
    ("pdb_id / best_chain", "PDB entry and the chain with the most observed residues in the aligned region."),
    ("identity_pct", "Identical residues / alignment columns (incl. gaps). Main predictor of model accuracy."),
    ("similarity_pct", "Residue pairs with positive BLOSUM62 score / alignment columns."),
    ("identity_vs_query_len_pct", "Identical residues / full query length (penalises partial hits)."),
    ("query_coverage_pct", "% of query residues aligned to a template residue. Residues not covered must be modelled de novo."),
    ("template_coverage_pct", "% of template sequence used by the alignment."),
    ("query_start/end, template_start/end", "Alignment boundaries (template numbering = PDB entity/SEQRES numbering)."),
    ("gap_opens / gaps_in_*", "Number of indels; many gaps = harder alignment, loop modelling needed."),
    ("evalue / bitscore", "Search statistics. E-value < 1e-3 is a confident homolog."),
    ("hssp_threshold_pct / hssp_distance", "Rost (1999) twilight-zone curve for this alignment length; distance > 0 = homology reliably inferred from identity."),
    ("method / resolution", "Experimental method and resolution (Å). Lower is better; ≤2.5 Å preferred."),
    ("r_free / r_work", "Crystallographic agreement factors; R-free ≤ 0.25 is good."),
    ("clashscore / rama_outliers_pct / rotamer_outliers_pct", "wwPDB validation geometry. Lower is better."),
    ("aligned_region_observed_pct", "% of the template residues used in the alignment that actually have coordinates (best chain)."),
    ("mutations", "Engineered mutations in the template entity."),
    ("oligomeric_state / stoichiometry", "Biological assembly; relevant if you need to model an oligomer."),
    ("ligands / metal_ions", "Bound ligands (common crystallisation additives removed) and metal ions; choose holo vs apo template accordingly."),
    ("warnings", "Automatic flags that deserve a manual look."),
    ("score_* columns", "0-1 component scores: identity (20→100%), similarity (30→100%), coverage (0→100%), resolution (3.5→1.5 Å; NMR=0.35), completeness (70→100% observed), r_free (0.32→0.20), geometry (clashscore 40→5, Rama outliers 5→0.5%). Missing data = 0.5."),
]


def write_excel(path, query_info, cfg, df, uniq, multi, coverage_sets, aln_store, uncovered, run_meta):
    from openpyxl import Workbook
    from openpyxl.formatting.rule import ColorScaleRule
    from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
    from openpyxl.utils import get_column_letter

    HDR = PatternFill("solid", fgColor="1F4E78")
    HFONT = Font(bold=True, color="FFFFFF")
    TITLE = Font(bold=True, size=14, color="1F4E78")
    BOLD = Font(bold=True)
    MONO = Font(name="Consolas", size=9)
    thin = Side(style="thin", color="D9D9D9")
    BORDER = Border(top=thin, bottom=thin, left=thin, right=thin)

    def green_good(ws, col_letter, n, low_good=False):
        rng = f"{col_letter}2:{col_letter}{n + 1}"
        if low_good:
            ws.conditional_formatting.add(rng, ColorScaleRule(start_type="min", start_color="63BE7B",
                                                              mid_type="percentile", mid_value=50, mid_color="FFEB84",
                                                              end_type="max", end_color="F8696B"))
        else:
            ws.conditional_formatting.add(rng, ColorScaleRule(start_type="min", start_color="F8696B",
                                                              mid_type="percentile", mid_value=50, mid_color="FFEB84",
                                                              end_type="max", end_color="63BE7B"))

    def write_df(ws, frame, start_row=1, widths=None):
        cols = list(frame.columns)
        for j, c in enumerate(cols, 1):
            cell = ws.cell(row=start_row, column=j, value=c)
            cell.fill, cell.font = HDR, HFONT
            cell.alignment = Alignment(wrap_text=True, vertical="center")
        for i, rec in enumerate(frame.itertuples(index=False), start_row + 1):
            for j, v in enumerate(rec, 1):
                if isinstance(v, (np.floating, float)) and (v != v):
                    v = None
                elif isinstance(v, np.generic):
                    v = v.item()
                ws.cell(row=i, column=j, value=v).border = BORDER
        for j, c in enumerate(cols, 1):
            L = get_column_letter(j)
            if widths and c in widths:
                ws.column_dimensions[L].width = widths[c]
            else:
                mx = max([len(str(c))] + [len(str(x)) for x in frame[c].head(200).tolist()])
                ws.column_dimensions[L].width = min(max(9, mx + 2), 45)
        ws.row_dimensions[start_row].height = 32
        return cols

    def link_pdb(ws, cols, n, start_row=1):
        if "pdb_id" not in cols:
            return
        j = cols.index("pdb_id") + 1
        for i in range(start_row + 1, start_row + n + 1):
            c = ws.cell(row=i, column=j)
            if c.value:
                c.hyperlink = f"https://www.rcsb.org/structure/{c.value}"
                c.font = Font(color="0563C1", underline="single")

    def style_table(ws, cols, n):
        ws.freeze_panes = ws.cell(row=2, column=min(4, len(cols)) + 1)
        ws.auto_filter.ref = f"A1:{get_column_letter(len(cols))}{n + 1}"
        good = ["template_score", "identity_pct", "similarity_pct", "query_coverage_pct",
                "aligned_region_observed_pct", "bitscore", "hssp_distance"]
        bad = ["resolution", "r_free", "clashscore", "rama_outliers_pct"]
        for c in good + bad:
            if c in cols and n > 1:
                green_good(ws, get_column_letter(cols.index(c) + 1), n, low_good=c in bad)
        for c, fmt in [("evalue", "0.00E+00"), ("resolution", "0.00"), ("r_free", "0.000"),
                       ("r_work", "0.000"), ("bitscore", "0.0")]:
            if c in cols:
                L = get_column_letter(cols.index(c) + 1)
                for i in range(2, n + 2):
                    ws[f"{L}{i}"].number_format = fmt
        if "template_score" in cols:
            L = get_column_letter(cols.index("template_score") + 1)
            for i in range(2, n + 2):
                ws[f"{L}{i}"].font = BOLD

    wb = Workbook()

    # ---- Summary
    ws = wb.active
    ws.title = "Summary"
    ws["A1"] = "Homology modelling template report"
    ws["A1"].font = TITLE
    info_rows = [
        ("Query", query_info["name"]),
        ("Source", query_info["source"]),
        ("Query length", len(query_info["sequence"])),
        ("Search engine", cfg["search_engine"].upper()),
        ("E-value cutoff", cfg["evalue_cutoff"]),
        ("Hits (polymer entities)", len(df)),
        ("Distinct PDB entries", df.pdb_id.nunique() if len(df) else 0),
        ("Unique template sequences", len(uniq)),
        ("Hits in tier A / B / C / D",
         " / ".join(str(int((df.confidence_tier.str[0] == t).sum())) if len(df) else "0" for t in "ABCD")),
        ("Query residues covered by ≥1 hit",
         f"{len(set().union(*coverage_sets.values())) if coverage_sets else 0} / {len(query_info['sequence'])}"),
        ("Regions with NO template",
         ", ".join(f"{a}-{b}" for a, b in uncovered) if uncovered else "none (full coverage)"),
        ("Run date", run_meta["date"]),
        ("Run time (s)", run_meta["seconds"]),
    ]
    for i, (k, v) in enumerate(info_rows, 3):
        ws.cell(row=i, column=1, value=k).font = BOLD
        ws.cell(row=i, column=2, value=v)
    r0 = 3 + len(info_rows) + 1
    ws.cell(row=r0, column=1, value="Query sequence").font = BOLD
    seq = query_info["sequence"]
    for k in range(0, len(seq), 60):
        c = ws.cell(row=r0 + k // 60, column=2, value=f"{k + 1:>5}  {seq[k:k + 60]}")
        c.font = MONO
    r = r0 + (len(seq) - 1) // 60 + 2

    ws.cell(row=r, column=1, value="Top 10 recommended templates (unique sequences)").font = TITLE
    if len(uniq):
        top = uniq.head(10)[["unique_rank", "pdb_id", "best_chain", "template_score", "confidence_tier",
                             "identity_pct", "query_coverage_pct", "query_start", "query_end", "resolution",
                             "method", "organism", "description", "ligands", "warnings"]]
        cols = write_df(ws, top, start_row=r + 1)
        link_pdb(ws, cols, len(top), start_row=r + 1)
        r = r + len(top) + 3
    if len(multi) > 1:
        ws.cell(row=r, column=1, value="Suggested multi-template set (complementary coverage)").font = TITLE
        write_df(ws, multi, start_row=r + 1)
        r += len(multi) + 3
    ws.cell(row=r, column=1, value="How to choose").font = TITLE
    tips = [
        "1. Start from the top of 'Unique_Templates'; prefer tier A/B, coverage of your region of interest, and good resolution.",
        "2. Check 'Coverage_Map' – if no single template covers the whole query, use the multi-template set or model domains separately.",
        "3. Match the biological state you need: ligand-bound (holo) vs apo, oligomeric state, active/inactive conformation.",
        "4. Avoid templates with many missing residues, mutations at key sites, or 'below HSSP curve' warnings.",
        "5. Inspect the alignment ('Alignments' sheet) around gaps and functional residues before modelling.",
    ]
    for i, t in enumerate(tips, r + 1):
        ws.cell(row=i, column=1, value=t)
    ws.column_dimensions["A"].width = max(ws.column_dimensions["A"].width or 0, 34)
    ws.column_dimensions["B"].width = max(ws.column_dimensions["B"].width or 0, 75)

    # ---- Ranked hits & unique templates
    order = ["rank", "pdb_id", "best_chain", "entity_id", "template_score", "confidence_tier", "identity_pct",
             "similarity_pct", "query_coverage_pct", "identity_vs_query_len_pct", "query_start", "query_end",
             "template_start", "template_end", "template_length", "template_coverage_pct", "aligned_residues",
             "gap_opens", "gaps_in_query", "gaps_in_template", "evalue", "bitscore", "hssp_threshold_pct",
             "hssp_distance", "homology_reliable", "method", "resolution", "r_free", "r_work", "clashscore",
             "rama_outliers_pct", "rotamer_outliers_pct", "aligned_region_observed_pct", "mutations",
             "oligomeric_state", "stoichiometry", "n_protein_entities", "n_chains_entry", "ligands", "metal_ions",
             "description", "organism", "uniprot_ids", "same_organism", "release_date", "title", "warnings",
             "all_chains"] + [c for c in df.columns if c.startswith("score_")]
    if not cfg.get("prefer_organism"):
        order.remove("same_organism")
    for name, frame in [("Ranked_Hits", df), ("Unique_Templates", uniq)]:
        ws = wb.create_sheet(name)
        if frame.empty:
            ws["A1"] = "No hits."
            continue
        cols_order = (["unique_rank"] if "unique_rank" in frame.columns else []) + order + (
            ["n_structures_same_sequence", "other_pdb_ids_same_sequence"]
            if "n_structures_same_sequence" in frame.columns else [])
        out = frame[[c for c in cols_order if c in frame.columns]]
        cols = write_df(ws, out, widths={"title": 50, "description": 35, "warnings": 45, "ligands": 35,
                                          "confidence_tier": 30, "other_pdb_ids_same_sequence": 40})
        link_pdb(ws, cols, len(out))
        style_table(ws, cols, len(out))

    # ---- Coverage map
    ws = wb.create_sheet("Coverage_Map")
    qlen = len(query_info["sequence"])
    b = cfg["coverage_bin"]
    nb = math.ceil(qlen / b)
    ws["A1"] = (f"Fraction of each {b}-residue window of the query covered by each template "
                "(green = covered, red = not covered). Top unique templates, by rank.")
    ws["A1"].font = BOLD
    hdr = ["rank", "pdb_id", "chain", "identity_pct", "coverage_pct"] + [f"{i * b + 1}-{min((i + 1) * b, qlen)}"
                                                                         for i in range(nb)]
    for j, h in enumerate(hdr, 1):
        c = ws.cell(row=3, column=j, value=h)
        c.fill, c.font = HDR, HFONT
        c.alignment = Alignment(text_rotation=90 if j > 5 else 0, horizontal="center")
    green, red = PatternFill("solid", fgColor="63BE7B"), PatternFill("solid", fgColor="F8696B")
    yellow = PatternFill("solid", fgColor="FFEB84")
    rows_cm = list(uniq.head(cfg["n_coverage_map"]).itertuples()) if len(uniq) else []
    union = set()
    for i, rr in enumerate(rows_cm, 4):
        cs = coverage_sets[rr.entity_id]
        union |= cs
        vals = [rr.unique_rank, rr.pdb_id, rr.best_chain, rr.identity_pct, rr.query_coverage_pct]
        for j, v in enumerate(vals, 1):
            ws.cell(row=i, column=j, value=v)
        ws.cell(row=i, column=2).hyperlink = f"https://www.rcsb.org/structure/{rr.pdb_id}"
        for k in range(nb):
            win = range(k * b + 1, min((k + 1) * b, qlen) + 1)
            frac = sum(1 for p in win if p in cs) / len(win)
            c = ws.cell(row=i, column=6 + k, value=round(frac, 2))
            c.fill = green if frac >= 0.8 else (yellow if frac > 0 else red)
            c.font = Font(size=7)
    rU = 4 + len(rows_cm)
    ws.cell(row=rU, column=2, value="UNION (all hits)").font = BOLD
    allcov = set().union(*coverage_sets.values()) if coverage_sets else set()
    for k in range(nb):
        win = range(k * b + 1, min((k + 1) * b, qlen) + 1)
        frac = sum(1 for p in win if p in allcov) / len(win)
        c = ws.cell(row=rU, column=6 + k, value=round(frac, 2))
        c.fill = green if frac >= 0.8 else (yellow if frac > 0 else red)
        c.font = Font(size=7, bold=True)
    for j in range(6, 6 + nb):
        ws.column_dimensions[get_column_letter(j)].width = 4
    ws.column_dimensions["B"].width = 16
    ws.freeze_panes = "F4"

    # ---- Alignments
    ws = wb.create_sheet("Alignments")
    ws["A1"] = "Pairwise alignments (| identical, : similar). Template numbering = PDB SEQRES/entity numbering."
    ws["A1"].font = BOLD
    r = 3
    width = 60
    for rr in uniq.head(cfg["n_alignments_in_excel"]).itertuples() if len(uniq) else []:
        qa, sa = aln_store[rr.entity_id]
        ml = match_line(qa, sa)
        ws.cell(row=r, column=1, value=f"#{rr.unique_rank}  {rr.pdb_id} chain {rr.best_chain}  |  "
                                       f"id {rr.identity_pct}%  cov {rr.query_coverage_pct}%  "
                                       f"E={(rr.evalue if rr.evalue is not None and rr.evalue == rr.evalue else float('nan')):.2g}  |  {rr.description}").font = BOLD
        r += 1
        qi, si = rr.query_start, rr.template_start
        for k in range(0, len(qa), width):
            qs, ss, ms = qa[k:k + width], sa[k:k + width], ml[k:k + width]
            ws.cell(row=r, column=1, value=f"Query  {qi:>5} {qs}").font = MONO
            ws.cell(row=r + 1, column=1, value=f"             {ms}").font = MONO
            ws.cell(row=r + 2, column=1, value=f"{rr.pdb_id}_{str(rr.best_chain):<2} {si:>5} {ss}").font = MONO
            qi += sum(1 for x in qs if x != "-")
            si += sum(1 for x in ss if x != "-")
            r += 4
        r += 1
    ws.column_dimensions["A"].width = 90

    # ---- FASTA-style alignment for modelling tools
    ws = wb.create_sheet("Alignment_FASTA")
    ws["A1"] = "Aligned FASTA pairs (copy into MODELLER/SWISS-MODEL 'user alignment' after checking)."
    ws["A1"].font = BOLD
    r = 3
    for rr in uniq.head(cfg["n_alignments_in_excel"]).itertuples() if len(uniq) else []:
        qa, sa = aln_store[rr.entity_id]
        ws.cell(row=r, column=1, value=f">query/{rr.query_start}-{rr.query_end}").font = MONO
        ws.cell(row=r + 1, column=1, value=qa).font = MONO
        ws.cell(row=r + 2, column=1, value=f">{rr.pdb_id}_{rr.best_chain}/{rr.template_start}-{rr.template_end}").font = MONO
        ws.cell(row=r + 3, column=1, value=sa).font = MONO
        r += 5
    ws.column_dimensions["A"].width = 120

    # ---- Metric guide & settings
    ws = wb.create_sheet("Metric_Guide")
    ws.append(["Column", "Meaning"])
    for c in ws[1]:
        c.fill, c.font = HDR, HFONT
    for k, v in COLUMN_DOCS:
        ws.append([k, v])
    for row in ws.iter_rows(min_row=2):
        row[1].alignment = Alignment(wrap_text=True, vertical="top")
        row[0].font = BOLD
    ws.column_dimensions["A"].width = 42
    ws.column_dimensions["B"].width = 110

    ws = wb.create_sheet("Settings")
    ws.append(["Parameter", "Value"])
    for c in ws[1]:
        c.fill, c.font = HDR, HFONT
    for k, v in cfg.items():
        if k == "weights":
            for wk, wv in v.items():
                ws.append([f"weight_{wk}", wv])
        else:
            ws.append([k, json.dumps(v) if isinstance(v, (list, dict)) else v])
    ws.column_dimensions["A"].width = 30
    ws.column_dimensions["B"].width = 40

    wb.save(path)


# ---------------------------------------------------------------------------
# 8. Structure download
# ---------------------------------------------------------------------------
def download_structures(pdb_ids, zip_path):
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
        for pid in pdb_ids:
            for ext in ("pdb", "cif"):
                r = _request("GET", f"https://files.rcsb.org/download/{pid}.{ext}")
                if r.status_code == 200 and r.content:
                    z.writestr(f"{pid}.{ext}", r.content)
                    break
    return zip_path


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def run_pipeline(**overrides):
    cfg = json.loads(json.dumps(DEFAULT_CONFIG))
    for k, v in overrides.items():
        if k == "weights" and v:
            cfg["weights"].update(v)
        elif v is not None:
            cfg[k] = v
    t0 = time.time()

    print("[1/5] Resolving query ...")
    qi = fetch_query_sequence(cfg["query"])
    qseq = qi["sequence"]
    if len(qseq) < 20:
        raise ValueError("Query sequence must be at least 20 residues.")
    print(f"  {qi['source']}: {qi['name'][:80]}  ({len(qseq)} aa)")

    print(f"[2/5] Searching PDB for homologs ({cfg['search_engine']}) ...")
    if cfg["search_engine"].lower() == "ncbi":
        hits = search_ncbi(qseq, cfg["evalue_cutoff"], cfg["max_hits"])
    else:
        hits = search_rcsb(qseq, cfg["evalue_cutoff"], cfg["min_identity_pct"] / 100.0, cfg["max_hits"])
    print(f"  {len(hits)} homologous polymer entities found")
    if not hits:
        print("No homologs found. Try a looser E-value cutoff or the NCBI search engine.")
        return None

    print("[3/5] Fetching structure metadata & validation metrics ...")
    meta = fetch_entity_metadata([h["entity_id"] for h in hits])

    print("[4/5] Computing alignment metrics and ranking ...")
    df, cov_sets, aln_store = build_table(qi, hits, meta, cfg)
    if df.empty:
        print("All hits were removed by the filters.")
        return None
    cov_sets = {k: v for k, v in cov_sets.items() if k in set(df.entity_id)}
    uniq = unique_templates(df)
    multi = greedy_multitemplate(uniq, cov_sets, len(qseq))
    unc = uncovered_regions(list(cov_sets.values()), len(qseq))

    print("[5/5] Writing Excel ...")
    safe_name = re.sub(r"[^A-Za-z0-9]+", "_", qi["source"].split()[-1] if qi["source"] != "User sequence"
                       else "sequence")[:30]
    stamp = datetime.now().strftime("%Y%m%d_%H%M")
    xlsx = f"{cfg['output_prefix']}_{safe_name}_{stamp}.xlsx"
    write_excel(xlsx, qi, cfg, df.drop(columns=["_entity_sequence"]), uniq.drop(columns=["_entity_sequence"]),
                multi, cov_sets, aln_store, unc,
                {"date": datetime.now().strftime("%Y-%m-%d %H:%M"), "seconds": round(time.time() - t0, 1)})
    outputs = {"excel": xlsx}

    if cfg["download_top_n"] and len(uniq):
        zp = f"{cfg['output_prefix']}_{safe_name}_top{cfg['download_top_n']}_structures.zip"
        print(f"  Downloading top {cfg['download_top_n']} template structures ...")
        outputs["structures_zip"] = download_structures(list(dict.fromkeys(uniq.pdb_id))[:cfg["download_top_n"]], zp)

    print(f"\nDone in {time.time() - t0:.0f}s. {len(df)} hits, {len(uniq)} unique templates.")
    cols = ["unique_rank", "pdb_id", "best_chain", "template_score", "identity_pct", "query_coverage_pct",
            "resolution", "method", "description"]
    print(uniq[cols].head(10).to_string(index=False))
    if unc:
        print("\nQuery regions with no template:", ", ".join(f"{a}-{b}" for a, b in unc))
    return {"outputs": outputs, "ranked": df, "unique": uniq, "multi_template": multi, "query": qi}


if __name__ == "__main__":
    import sys
    run_pipeline(query=sys.argv[1] if len(sys.argv) > 1 else DEFAULT_CONFIG["query"])
