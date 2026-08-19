"""
PRKN Follow-up: Isolate the true catalytic residue (Cys431) from the
structural zinc-ligand scaffold, to test whether the depletion signal
is driven specifically by the zinc-finger residues (as hypothesized)
rather than by the catalytic site itself.
"""

import pandas as pd
import numpy as np
from scipy import stats

PROJECT_DIR = 'pathogenicity_abundance_project'

# Load the already-computed PRKN results (from the main pipeline run)
valid = pd.read_csv(f'{PROJECT_DIR}/PRKN_final_results.csv')

# Re-derive is_failure_case if not already present (standard thresholds)
if 'is_failure_case' not in valid.columns:
    valid['is_failure_case'] = (valid['a2a_prediction'] < 0.4) & (valid['score'] > 0.8)

def extract_position(variant):
    import re
    m = re.match(r'^[A-Z](\d+)[A-Z]$', str(variant))
    return int(m.group(1)) if m else None

if 'position' not in valid.columns:
    valid['position'] = valid['variant_1letter'].apply(extract_position)

# ============================================================
# Define three site sets to compare
# ============================================================

SITE_SETS = {
    'full_override_25pos': [431, 238, 241, 253, 257, 260, 263, 289, 293, 332, 337,
                              352, 360, 365, 368, 373, 377, 418, 421, 436, 441,
                              446, 449, 457, 461],
    'catalytic_only_C431': [431],
    'ring2_region_excl_zinc': [431, 418, 421, 436, 441, 446, 449, 457, 461],
    'confirmed_ibr_zinc_ligands_only': [332, 337, 352, 360, 365, 368, 373, 377],
}

print(f"{'Site set':<30} {'N sites':>8} {'N variants near':>16} "
      f"{'OR':>8} {'p-value':>12} {'Contingency table'}")
print("-" * 100)

for set_name, positions in SITE_SETS.items():
    valid['near_site_test'] = valid['position'].isin(positions)
    n_near = valid['near_site_test'].sum()

    contingency = pd.crosstab(valid['is_failure_case'], valid['near_site_test'])

    if contingency.shape == (2, 2) and n_near >= 3:
        odds_ratio, p_value = stats.fisher_exact(contingency)
    else:
        odds_ratio, p_value = np.nan, np.nan

    print(f"{set_name:<30} {len(positions):>8} {n_near:>16} "
          f"{odds_ratio:>8.2f} {p_value:>12.2e}")
    print(f"  Contingency table:\n{contingency}\n")

print("\nINTERPRETATION GUIDE:")
print("- If 'catalytic_only_C431' shows OR > 1 (even if not significant, n=1 position)")
print("  while 'confirmed_ibr_zinc_ligands_only' shows OR << 1, this confirms the")
print("  zinc-ligand hypothesis: the depletion is driven by structural residues,")
print("  not the catalytic site itself.")
print("- If both show depletion, the explanation needs revisiting.")
