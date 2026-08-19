# A2A: relating AlphaMissense pathogenicity predictions to experimental protein abundance

Analysis code and derived data for the manuscript:

> Azza Althagafi. *Evidence for functional-site decoupling in protein abundance
> assays using AlphaMissense scores.* Manuscript in preparation.

A single fixed transformation of AlphaMissense pathogenicity scores is compared
against multiplexed measurements of protein abundance and stability for six
proteins (PTEN, TPMT, NUDT15, ASPA, PRKN, VKORC1), and variants predicted
pathogenic yet experimentally stable are tested for positional enrichment at
annotated catalytic and functional sites.

---

## Quick start

```bash
git clone https://github.com/TODO_USER/a2a-abundance.git
cd a2a-abundance
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

python code/download_data.py         # fetch MaveDB scores + AlphaMissense extracts
python code/unified_pipeline.py      # main analysis, writes results/
python code/prkn_site_isolation.py   # PRKN catalytic vs. structural decomposition
python code/make_all_figures.py      # figures 1-2 and supplementary figure S1
python code/verify_prkn.py results/PRKN_final_results.csv   # reproduces Table 4
```

The committed inputs under `data/` are sufficient to reproduce every number in
the manuscript without any network access. See "Reproducibility" below.

---

## Repository contents

| Path | Contents |
|---|---|
| `code/download_data.py` | Downloads MaveDB score sets and extracts per-protein AlphaMissense predictions; documents which score set was chosen for each protein and why |
| `code/unified_pipeline.py` | Main six-protein analysis: matching, A2A transform, correlation, enrichment testing, threshold sensitivity, FDR correction |
| `code/prkn_site_isolation.py` | Decomposition of PRKN functional sites into catalytic and structural subsets |
| `code/make_all_figures.py` | All manuscript figures |
| `code/verify_prkn.py` | Standalone check reproducing the PRKN decomposition (manuscript Table 4) |
| `data/mavedb/` | Experimental score sets as downloaded from MaveDB |
| `data/alphamissense/` | Per-protein AlphaMissense extracts |
| `results/` | Per-variant results, threshold grids, summary table, Supplementary Table S1 |
| `figures/` | Generated figures (PDF and PNG) |

---

## Data sources

**MaveDB** — experimental abundance and stability scores. The score sets used are:

| Protein | MaveDB accession |
|---|---|
| PTEN | `urn:mavedb:00000013-a-1` |
| TPMT | `urn:mavedb:00000013-b-1` |
| NUDT15 | `urn:mavedb:00000055-a-1` |
| ASPA | `urn:mavedb:00000657-a-1` |
| PRKN | `urn:mavedb:00000114-a-1` |
| VKORC1 | `urn:mavedb:00000078-b-1` |

**AlphaMissense** — predictions from `AlphaMissense_aa_substitutions.tsv.gz`,
available from https://github.com/google-deepmind/alphamissense. This file is
~1 GB and is **not** included here. It is only needed to regenerate the
per-protein extracts in `data/alphamissense/`, which are committed; place it in
`data/alphamissense/` if you wish to re-extract from source.

**UniProt** — functional-site annotations retrieved via the REST API,
**release 2026_02 (10 June 2026)**. Because UniProt annotations change between
releases, re-running the retrieval today may return a different set of features
than reported in the manuscript; `results/Supplementary_Table_S1.csv` records
the retrieval as used.

---

## Reproducibility

Two things are worth knowing before re-running.

**Data provenance.** `download_data.py` fetches everything from source;
`unified_pipeline.py` will also download on demand if the caches are absent.
Both write to the same locations, so either entry point works.


**The pipeline queries live APIs.** `unified_pipeline.py` downloads from MaveDB
and UniProt on first run and caches the results. The cached files are committed,
so a fresh clone reproduces the published numbers exactly. Deleting them will
cause the script to re-download, and results may differ if either resource has
been updated.

**UniProt drift is expected.** The manuscript documents a case where automated
UniProt annotation spans entire structural domains rather than discrete
catalytic residues, and where curation reverses the direction of the apparent
effect (PRKN). The curated site definitions are hardcoded in
`MANUAL_SITE_OVERRIDES` in `unified_pipeline.py`, with the literature basis for
each documented in Supplementary Table S2 of the manuscript.


To confirm a clean reproduction, check that `results/all_proteins_summary.csv`
matches:

| Protein | Spearman ρ | Enrichment OR | p |
|---|---|---|---|
| PTEN | 0.482 | 2.64 | 1.19 × 10⁻¹⁵ |
| TPMT | 0.552 | 1.97 | 2.47 × 10⁻⁴ |
| NUDT15 | 0.487 | 2.19 | 2.18 × 10⁻⁸ |
| ASPA | 0.486 | 15.53 | 3.26 × 10⁻⁵⁴ |
| PRKN | 0.637 | 0.16 | 3.40 × 10⁻¹⁶ |
| VKORC1 | 0.483 | 4.11 | 4.34 × 10⁻⁸ |

---

## Citation

This manuscript is in preparation and not yet published; a full citation will
be added here on publication. In the meantime, please cite this repository
directly, and note the commit or release you used.

Please also cite the underlying resources: MaveDB (Esposito et al. 2019),
AlphaMissense (Cheng et al. 2023), and UniProt (The UniProt Consortium 2025).

## License

Code is released under the MIT License (see `LICENSE`). Derived data in
`results/` is released under CC BY 4.0. The original MaveDB, AlphaMissense, and
UniProt data remain subject to their respective licences.
