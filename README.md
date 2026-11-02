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
2. Searches the PDB for homologs: RCSB sequence search (MMseqs2, default) or NCBI BLASTp against `pdb`.
3. Fetches structure metadata from the RCSB GraphQL API: method, resolution, R-free/R-work, clashscore, Ramachandran and rotamer outliers, unobserved residues per chain, ligands (common additives filtered), metal ions, oligomeric state, mutations, organism, release date.
4. Computes alignment metrics: identity, BLOSUM62 similarity, query and template coverage, gaps, E-value, bit score, and distance above the Rost (1999) HSSP twilight-zone curve.
5. Ranks hits with a 0–100 composite template score (weights adjustable), assigns confidence tiers A–D, and flags issues.

## Excel output

| Sheet | Contents |
|---|---|
| Summary | Query info, top 10 templates, uncovered query regions, suggested multi-template set, selection tips |
| Ranked_Hits | Every hit with all metrics, colour scales, filters, RCSB links |
| Unique_Templates | Best structure per distinct template sequence |
| Coverage_Map | Per-window coverage of the query by each template |
| Alignments / Alignment_FASTA | Readable pairwise alignments and FASTA pairs for MODELLER / SWISS-MODEL |
| Metric_Guide / Settings | Column definitions and run parameters |

## Default score weights

| Component | Weight | Scaling (0 → 1) |
|---|---|---|
| Identity | 0.30 | 20% → 100% |
| Coverage | 0.25 | 0% → 100% |
| Resolution | 0.15 | 3.5 Å → 1.5 Å (NMR = 0.35) |
| Completeness | 0.10 | 70% → 100% of aligned template residues observed |
| R-free | 0.075 | 0.32 → 0.20 |
| Geometry | 0.075 | clashscore 40 → 5, Rama outliers 5% → 0.5% |
| Similarity | 0.05 | 30% → 100% |

Missing data scores 0.5.

## Requirements

`requests`, `pandas`, `numpy`, `openpyxl`; `biopython` for the NCBI search option (installed by the notebook).
