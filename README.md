# Homology Template Finder

Find every experimentally solved homolog of a protein in the PDB and rank them as templates for homology modelling. Runs in Google Colab and writes a formatted Excel report.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/victorolatunde23-dev/homology_modelling_task/blob/main/Homology_Template_Finder.ipynb)

## Usage

1. Open `Homology_Template_Finder.ipynb` in Colab (badge above, or File ▸ Upload notebook).
2. Enter a query in the **Inputs** cell: a UniProt accession (`P69905`), a PDB ID with optional chain (`4HHB_A`), an amino-acid sequence, or upload a FASTA file. Several queries can be comma-separated.
3. Runtime ▸ Run all. The Excel file (and optional structure zip) downloads automatically.

`homology_template_finder.py` is the same code as a standalone module:

```python
from homology_template_finder import run_pipeline
res = run_pipeline(query="P69905", evalue_cutoff=1e-3, download_top_n=5)
```

## What it does

1. Resolves the query sequence (UniProt REST or RCSB FASTA).
2. Searches the PDB for homologs with one of three engines:
   - `rcsb` (default): RCSB sequence search (MMseqs2, pairwise), fast.
   - `ncbi`: NCBI BLASTp against `pdb`.
   - `colabfold`: the ColabFold MMseqs2 server builds an MSA of the query (UniRef + environmental databases) and searches the PDB with that profile, the same template search ColabFold uses for AlphaFold2. It finds remote homologs a pairwise search misses and reports how many sequences the MSA contains. It takes a few minutes and runs on a shared server (~20–50k requests/day), so use it for individual queries rather than bulk runs.

   If fewer than 5 hits come back, or at least 25% of the query is low-complexity, the tool re-searches NCBI with the low-complexity filter and composition-based statistics switched off. Masking of this kind caused ColabFold's one notable CASP14 miss (T1084). Hits found only this way are marked `found_by = ncbi-unmasked`.
3. Fetches structure metadata from the RCSB GraphQL API: method, resolution, R-free/R-work, clashscore, Ramachandran and rotamer outliers, unobserved residues per chain, ligands (common additives filtered), metal ions, oligomeric state, mutations, organism, release date.
4. Computes alignment metrics: identity, BLOSUM62 similarity, query and template coverage, gaps, E-value, bit score, and distance above the Rost (1999) HSSP twilight-zone curve.
5. Ranks hits with a 0–100 composite template score (weights adjustable), assigns confidence tiers A–D, and flags issues (including homo-oligomeric templates).
6. Groups near-identical templates: a template ≥95% identical to a better-ranked one, over the query positions both cover, joins its cluster. The Summary lists cluster representatives only, plus the best templates in each identity-to-query band (0–20, 20–40, … %), so remote but structurally different templates are not crowded out. This follows ColabFold's diversity-aware MSA filter.
7. If there is no tier A/B template, or part of the query has no template, writes a query CSV for the [ColabFold AlphaFold2 notebook](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb). Downloaded templates are lower-case mmCIF files (`4hhb.cif`) that ColabFold accepts as custom templates.

## Excel output

| Sheet | Contents |
|---|---|
| Summary | Query info, notes (MSA depth, unmasked hits, ColabFold hand-off), top 10 non-redundant templates, best templates per identity band, uncovered query regions, suggested multi-template set, selection tips |
| Ranked_Hits | Every hit with all metrics, colour scales, filters, RCSB links |
| Unique_Templates | Best structure per distinct template sequence, with its redundancy cluster |
| Coverage_Map | Per-window coverage of the query by each template |
| Alignments / Alignment_FASTA | Readable pairwise alignments and FASTA pairs for MODELLER / SWISS-MODEL |
| Metric_Guide / Settings | Column definitions and run parameters |

## Default score weights

| Component | Weight | Scaling (0 → 1) |
|---|---|---|
| Identity | 0.25 | 20% → 100% |
| Coverage | 0.20 | 0% → 100% |
| Significance | 0.10 | E-value 10 → 1e-30 |
| Resolution | 0.15 | 3.5 Å → 1.5 Å (NMR = 0.35) |
| Completeness | 0.10 | 70% → 100% of aligned template residues observed |
| R-free | 0.075 | 0.32 → 0.20 |
| Geometry | 0.075 | clashscore 40 → 5, Rama outliers 5% → 0.5% |
| Similarity | 0.05 | 30% → 100% |

Missing data scores 0.5. Significance was added so that a short local hit with high identity but a weak E-value no longer outranks a long, confident homolog.

## Other options

| Option | Default | Meaning |
|---|---|---|
| `max_release_date` | `""` | `YYYY-MM-DD`: ignore templates released later (benchmarking) |
| `redundancy_identity_pct` | 95 | Clustering threshold; 100 switches it off |
| `diverse_per_band` | 2 | Templates listed per identity band on the Summary sheet |
| `low_complexity_retry` / `low_complexity_min_hits` | `True` / 5 | Unmasked NCBI re-search rule |
| `colabfold_handoff` | `"auto"` | `"auto"`, `"always"` or `"never"` write the ColabFold query CSV |

## Tests

`pytest -q tests/` runs the whole pipeline offline with every web service mocked (RCSB, ColabFold server, NCBI).

## Requirements

`requests`, `pandas`, `numpy`, `openpyxl`; `biopython` for the NCBI search and the unmasked re-search (installed by the notebook).

## Reference

Mirdita M, Schütze K, Moriwaki Y, Heo L, Ovchinnikov S, Steinegger M. ColabFold: making protein folding accessible to all. *Nature Methods* 19, 679–682 (2022).
