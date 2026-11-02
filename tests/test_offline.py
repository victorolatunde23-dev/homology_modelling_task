"""Offline end-to-end tests: every web service is replaced by a mock.

Run with:  pytest -q tests/
"""
import io
import random
import sys
import tarfile
from pathlib import Path

import pytest
from openpyxl import load_workbook

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import homology_template_finder as htf  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"
rng = random.Random(7)
QUERY = "".join(rng.choice(AA) for _ in range(120))


def mutate(seq, frac, seed):
    r = random.Random(seed)
    out = list(seq)
    for i in r.sample(range(len(seq)), int(frac * len(seq))):
        out[i] = r.choice([a for a in AA if a != out[i]])
    return "".join(out)


# entity_id -> (template sequence, query region (1-based, inclusive), resolution, oligomer)
TEMPLATES = {
    "4AAA_1": (QUERY, (1, 120), 1.6, "Monomer"),
    "4AAB_1": (mutate(QUERY, 0.01, 1), (1, 120), 2.4, "Monomer"),      # 99% id -> same cluster as 4AAA
    "3BBB_1": (mutate(QUERY[:100], 0.55, 2), (1, 100), 2.0, "Homo 2-mer"),
    "2CCC_1": (mutate(QUERY[60:], 0.80, 3), (61, 120), 2.8, "Monomer"),
}


def fake_entity(eid, seq, res, olig):
    pdb = eid.split("_")[0]
    return {
        "rcsb_id": eid,
        "rcsb_polymer_entity": {"pdbx_description": f"protein {pdb}"},
        "entity_poly": {"pdbx_seq_one_letter_code_can": seq, "rcsb_sample_sequence_length": len(seq),
                        "rcsb_mutation_count": 0},
        "rcsb_entity_source_organism": [{"scientific_name": "Homo sapiens", "ncbi_taxonomy_id": 9606}],
        "rcsb_polymer_entity_container_identifiers": {"auth_asym_ids": ["A"], "uniprot_ids": []},
        "polymer_entity_instances": [{
            "rcsb_polymer_entity_instance_container_identifiers": {"auth_asym_id": "A", "asym_id": "A"},
            "rcsb_polymer_instance_feature": []}],
        "entry": {
            "rcsb_id": pdb, "struct": {"title": f"structure {pdb}"}, "exptl": [{"method": "X-RAY DIFFRACTION"}],
            "rcsb_entry_info": {"resolution_combined": [res], "polymer_entity_count_protein": 1,
                                "deposited_polymer_entity_instance_count": 2 if "2-mer" in olig else 1},
            "refine": [{"ls_R_factor_R_free": 0.22, "ls_R_factor_R_work": 0.18}],
            "rcsb_accession_info": {"deposit_date": "2015-01-01", "initial_release_date":
                                    "2024-06-01" if pdb == "2CCC" else "2016-01-01"},
            "nonpolymer_entities": [],
            "pdbx_vrpt_summary_geometry": [{"clashscore": 4, "percent_ramachandran_outliers": 0.2,
                                            "percent_rotamer_outliers": 1}],
            "assemblies": [{"rcsb_struct_symmetry": [{"kind": "Global Symmetry", "oligomeric_state": olig,
                                                      "stoichiometry": ["A2"] if "2-mer" in olig else ["A1"]}]}],
        },
    }


def fake_graphql(query, variables):
    ids = variables["ids"]
    if "entries(" in query:  # chain -> entity mapping
        return {"entries": [{"rcsb_id": i, "polymer_entities": [
            {"rcsb_id": f"{i}_1", "rcsb_polymer_entity_container_identifiers":
             {"auth_asym_ids": ["A"], "asym_ids": ["A"]}}]} for i in ids]}
    return {"polymer_entities": [fake_entity(e, TEMPLATES[e][0], TEMPLATES[e][2], TEMPLATES[e][3])
                                 for e in ids if e in TEMPLATES]}


class Resp:
    def __init__(self, status=200, js=None, content=b"", text=""):
        self.status_code, self._js, self.content, self.text = status, js, content, text

    def json(self):
        return self._js


def rcsb_result():
    out = []
    for eid, (seq, (qb, qe), _, _) in TEMPLATES.items():
        qa = QUERY[qb - 1:qe]
        sa = seq[:len(qa)]
        ident = sum(a == b for a, b in zip(qa, sa)) / len(qa)
        out.append({"identifier": eid, "services": [{"nodes": [{"match_context": [{
            "evalue": 10 ** (-60 * ident), "bitscore": 400 * ident, "sequence_identity": ident,
            "query_beg": qb, "query_end": qe, "subject_beg": 1, "subject_end": len(qa),
            "subject_length": len(seq), "query_aligned_seq": qa, "subject_aligned_seq": sa}]}]}]})
    return {"result_set": out}


def colabfold_tar():
    lines = []
    for eid, (seq, (qb, qe), _, _) in TEMPLATES.items():
        n = qe - qb + 1
        lines.append(f"101\t{eid.split('_')[0].lower()}_A\t0.5\t{n}\t0\t0\t{qb}\t{qe}\t1\t{n}\t1e-20\t150")
    files = {"pdb70.m8": "\n".join(lines) + "\n",
             "uniref.a3m": f">101\n{QUERY}\n" + "".join(f">s{i}\n{QUERY}\n" for i in range(12)),
             "bfd.mgnify30.metaeuk30.smag30.a3m": f">101\n{QUERY}\n>e1\n{QUERY}\n"}
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        for name, text in files.items():
            data = text.encode()
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    return buf.getvalue()


@pytest.fixture
def offline(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    calls = {"ncbi_unmasked": 0, "colabfold": []}

    def fake_request(method, url, **kw):
        if "search.rcsb.org" in url:
            return Resp(200, rcsb_result())
        if "files.rcsb.org" in url:
            return Resp(200, content=b"data_" + url.rsplit("/", 1)[1].encode() + b"\n_entity_poly_seq.entity_id 1\n")
        raise AssertionError(f"unexpected URL {url}")

    def fake_requests(method, url, **kw):
        calls["colabfold"].append(url)
        if url.endswith("ticket/msa"):
            assert kw["data"]["q"].startswith(">101\n")
            return Resp(200, {"status": "PENDING", "id": "abc"})
        if url.endswith("ticket/abc"):
            return Resp(200, {"status": "COMPLETE", "id": "abc"})
        if url.endswith("result/download/abc"):
            return Resp(200, content=colabfold_tar())
        raise AssertionError(url)

    def fake_ncbi(seq, evalue, max_hits, unmasked=False):
        calls["ncbi_unmasked"] += unmasked
        return [{"entity_id": "9ZZZ_1", "pdb_id": "9ZZZ"}] if unmasked else []

    monkeypatch.setattr(htf, "_request", fake_request)
    monkeypatch.setattr(htf, "_graphql", fake_graphql)
    monkeypatch.setattr(htf.requests, "request", fake_requests)
    monkeypatch.setattr(htf, "search_ncbi", fake_ncbi)
    monkeypatch.setattr(htf.time, "sleep", lambda s: None)
    return calls


def run(**kw):
    return htf.run_pipeline(query=QUERY, **kw)


def test_rcsb_ranking_clusters_and_excel(offline):
    res = run(download_top_n=3, low_complexity_min_hits=3)
    uniq = res["unique"]
    assert list(uniq.pdb_id[:2]) == ["4AAA", "4AAB"]
    # 99% identical template joins the cluster of the identical one
    assert uniq.set_index("pdb_id").loc["4AAB", "cluster_rep"] == uniq.set_index("pdb_id").loc["4AAA", "unique_rank"]
    reps = uniq[uniq.is_representative]
    assert set(reps.pdb_id) == {"4AAA", "3BBB", "2CCC"}
    # significance is part of the score
    assert "score_significance" in res["ranked"].columns
    # identity bands surface the remote template
    assert dict(zip(res["diverse"].pdb_id, res["diverse"].identity_band)) == {
        "4AAA": "80-100%", "3BBB": "40-60%", "2CCC": "20-40%"}
    # homo-oligomer flagged
    assert "homo-oligomer" in res["ranked"].set_index("pdb_id").loc["3BBB", "warnings"]
    # Excel has the new Summary sections and columns
    wb = load_workbook(res["outputs"]["excel"])
    summary = " ".join(str(c.value) for row in wb["Summary"].iter_rows() for c in row if c.value)
    assert "non-redundant" in summary and "Diverse templates by identity" in summary
    hdr = [c.value for c in wb["Unique_Templates"][1]]
    assert {"found_by", "cluster_rep", "n_in_cluster", "score_significance"} <= set(hdr)
    # structures: lower-case mmCIF, representatives only (4AAB is redundant)
    import zipfile
    names = zipfile.ZipFile(res["outputs"]["structures_zip"]).namelist()
    assert names == ["4aaa.cif", "3bbb.cif", "2ccc.cif"]
    # tier A template exists and the query is fully covered -> no ColabFold hand-off file
    assert "colabfold_input" not in res["outputs"]
    assert offline["ncbi_unmasked"] == 0


def test_max_release_date_and_handoff(offline):
    res = run(max_release_date="2020-01-01", max_identity_pct=60, low_complexity_min_hits=0)
    assert "2CCC" not in set(res["ranked"].pdb_id)          # released 2024
    assert set(res["ranked"].pdb_id) == {"3BBB"}
    # 3BBB covers 1-100 only, so 101-120 is uncovered -> ColabFold input written
    csv = Path(res["outputs"]["colabfold_input"]).read_text()
    assert csv.startswith("id,sequence\n") and QUERY in csv
    assert any("ColabFold" in n for n in res["notes"])


def test_low_complexity_retry(offline, monkeypatch):
    monkeypatch.setattr(htf, "search_rcsb", lambda *a, **k: [])
    res = htf.run_pipeline(query=QUERY)
    assert res is None                         # unmasked hit has no metadata in the mock
    assert offline["ncbi_unmasked"] == 1


def test_colabfold_engine(offline):
    res = run(search_engine="colabfold", low_complexity_min_hits=0)
    assert res["msa_info"]["msa_sequences"] == 13
    assert set(res["ranked"].found_by) == {"colabfold"}
    assert set(res["ranked"].pdb_id) == set(e.split("_")[0] for e in TEMPLATES)
    assert any("Only 13 homologous sequences" in n for n in res["notes"])
    assert any(u.endswith("ticket/msa") for u in offline["colabfold"])


HBA = ("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAH"
       "KLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR")


def test_low_complexity_fraction():
    assert htf.low_complexity_fraction(HBA) == 0.0                      # globular protein
    assert htf.low_complexity_fraction("MATLEKLMKAFESLKSF" + "Q" * 40 + "P" * 11) > 0.5   # polyQ/polyP


def test_redundancy_switch_off(offline):
    res = run(redundancy_identity_pct=100, low_complexity_min_hits=0)
    assert res["unique"].is_representative.all()
