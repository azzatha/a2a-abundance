#!/usr/bin/env python3
"""
Verify the PRKN decomposition numbers reported in Section 3.5 / Table 4.

Usage:
    python verify_prkn.py pathogenicity_abundance_project/PRKN_final_results.csv
"""

import sys
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact

CSV = sys.argv[1] if len(sys.argv) > 1 else "pathogenicity_abundance_project/PRKN_final_results.csv"

# Failure-case definition used throughout the manuscript
PRED_T, EXP_T = 0.4, 0.8

# PRKN functional-site subsets
CATALYTIC = [431]                                             # RING2 catalytic Cys
IBR       = [332, 337, 352, 360, 365, 368, 373, 377]          # IBR domain Zn ligands
RING0     = [238, 241, 253, 257, 260, 263, 289, 293]
RING2     = [418, 421, 436, 441, 446, 449, 457, 461]
ALL_ZN    = RING0 + IBR + RING2                               # 24 Zn-coordinating
AGGREGATE = CATALYTIC + ALL_ZN                                # 25 curated positions

df = pd.read_csv(CSV)
print(f"Loaded {len(df):,} rows from {CSV}\n")

# ---------------------------------------------------------------
# 1. Sanity checks against the file's own columns
# ---------------------------------------------------------------
print("=" * 62)
print("SANITY CHECKS")
print("=" * 62)

rule = (df.a2a_prediction < PRED_T) & (df.score > EXP_T)
match = bool((rule == df.is_failure_case).all())
print(f"is_failure_case == (a2a < {PRED_T} and score > {EXP_T}):  {match}"
      f"   [{int((rule != df.is_failure_case).sum())} mismatches]")

file_sites = set(df.loc[df.near_functional_site, "position"].unique())
print(f"near_functional_site positions: {len(file_sites)}  "
      f"(expected 25)   identical to hardcoded set: "
      f"{file_sites == set(AGGREGATE)}")
if file_sites != set(AGGREGATE):
    print(f"   only in file: {sorted(file_sites - set(AGGREGATE))}")
    print(f"   only in list: {sorted(set(AGGREGATE) - file_sites)}")

print(f"near_raw_uniprot_site positions: "
      f"{df.loc[df.near_raw_uniprot_site, 'position'].nunique()}  (expected 389)")
print(f"Overall failure-case rate: {df.is_failure_case.mean() * 100:.1f}%")

# ---------------------------------------------------------------
# 2. Enrichment per subset
# ---------------------------------------------------------------
def enrich(positions, label):
    """Fisher's exact test for one site subset vs the rest of the protein."""
    site = df.position.isin(positions)
    ct = (pd.crosstab(df.is_failure_case, site)
            .reindex(index=[False, True], columns=[False, True], fill_value=0))
    table = ct.values
    orr, p = fisher_exact(table)

    # Haldane-Anscombe correction when a cell is zero (OR/CI otherwise undefined)
    t = table + (0.5 if (table == 0).any() else 0.0)
    or_c = (t[1, 1] * t[0, 0]) / (t[1, 0] * t[0, 1])
    se = np.sqrt(sum(1.0 / x for x in t.flatten()))
    lo, hi = np.exp(np.log(or_c) - 1.96 * se), np.exp(np.log(or_c) + 1.96 * se)
    if (table == 0).any():
        lo = 0.0                      # zero observed failures -> lower bound 0

    n_var = int(site.sum())
    n_fail = int(df.loc[site, "is_failure_case"].sum())
    print(f"{label:<32}{len(positions):>5}{n_var:>7}{n_fail:>7}"
          f"{orr:>8.2f}  [{lo:.2f}-{hi:.2f}]  {p:.2e}")
    return orr, p

print("\n" + "=" * 62)
print("DECOMPOSITION (manuscript Table 4)")
print("=" * 62)
print(f"{'Subset':<32}{'pos':>5}{'vars':>7}{'fail':>7}{'OR':>8}  {'95% CI':<14}{'p':>9}")
enrich(AGGREGATE, "Full aggregate override")
enrich(CATALYTIC, "Catalytic (Cys431 only)")
enrich(IBR,       "Structural (IBR Zn ligands)")
enrich(ALL_ZN,    "Structural (all 24 Zn ligands)")
print("-" * 62)
enrich(RING0,     "  [RING0 Zn ligands]")
enrich(RING2,     "  [RING2 Zn ligands]")

# ---------------------------------------------------------------
# 3. The key claim: where do the failure cases actually sit?
# ---------------------------------------------------------------
print("\n" + "=" * 62)
print("KEY CLAIM")
print("=" * 62)
at_site = df[df.near_functional_site]
fails = at_site[at_site.is_failure_case]
print(f"Failure cases among the {len(at_site)} variants at annotated sites: {len(fails)}")
print(f"Of these, at Cys431: {int((fails.position == 431).sum())}")
print(f"At any Zn-coordinating residue: {int(fails.position.isin(ALL_ZN).sum())}")

c431 = df[df.position == 431]
print(f"\nCys431: {int(c431.is_failure_case.sum())}/{len(c431)} substitutions are "
      f"failure cases ({c431.is_failure_case.mean() * 100:.1f}%)")
print(f"Protein-wide rate: {df.is_failure_case.mean() * 100:.1f}%")