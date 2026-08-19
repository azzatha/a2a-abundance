"""
Pipeline steps:
  1. Load MaveDB experimental scores       <- ground truth
  2. Load AlphaMissense predictions        <- A2A input
  3. Match variants, apply A2A transform   <- core method
  4. Correlate A2A vs experimental score   <- Section 3.1 result
  5. Query UniProt functional sites        <- supplementary, for step 6 only
  6. Test failure-case enrichment at sites <- Section 3.4 result
"""

import pandas as pd
import numpy as np
import re
import requests
from pathlib import Path
from scipy import stats

PROJECT_DIR = Path('pathogenicity_abundance_project')

# ============================================================
# STEP 0: PROTEIN CONFIG 
# ============================================================

PROTEIN_CONFIGS = {
    'PTEN': {
        'mavedb_urn': 'urn:mavedb:00000013-a-1',
        'uniprot_id': 'P60484',
        'score_column': 'score',           
        'hgvs_column': 'hgvs_pro',         
    },
    'TPMT': {
        'mavedb_urn': 'urn:mavedb:00000013-b-1',
        'uniprot_id': 'P51580',
        'score_column': 'score',           
        'hgvs_column': 'hgvs_pro',
    },
    'NUDT15': {
        'mavedb_urn': 'urn:mavedb:00000055-a-1',   
        'uniprot_id': 'Q9NV35',
        'score_column': 'score',   
        'hgvs_column': 'hgvs_pro', 
    },
    # ============================================================
    # TIER 2 - Positive controls (proteins already examined by Livesey & Marsh 2025).
    # ============================================================
    'ASPA': {
        'mavedb_urn': 'urn:mavedb:00000657-a-1',                                                      
        'uniprot_id': 'P45381',
        'score_column': 'score',
        'hgvs_column': 'hgvs_pro',
    },
    'PRKN': {
        'mavedb_urn': 'urn:mavedb:00000114-a-1',  
        'uniprot_id': 'O60260',
        'score_column': 'score',
        'hgvs_column': 'hgvs_pro',
    },
    'VKORC1': {
        'mavedb_urn': 'urn:mavedb:00000078-b-1',   
        'uniprot_id': 'Q9BQB6',
        'score_column': 'score',
        'hgvs_column': 'hgvs_pro',
    },
}

# Known catalytic/functional sites from literature, used ONLY as a
# fallback when UniProt's automated annotation is too sparse (like PTEN).
# Document any entry here explicitly in the paper's supplementary table.
MANUAL_SITE_OVERRIDES = {
    'PTEN': {
        'positions': list(range(88, 99)) + list(range(123, 131)) + list(range(160, 172)),
        'source': 'WPD-loop (88-98), P-loop/HCxxGxxR motif (123-130), TI-loop (160-171); '
                   'multiple structural papers Structural Mechanisms '
                   'of PTEN Regulation',
        'reason': 'UniProt Active site annotation captured only C124; broader catalytic '
                   'architecture not annotated as Region/Motif/Binding site.'
    },
    'NUDT15': {
        'positions': list(range(47, 70)),  # Nudix box motif (48-69) + adjacent binding
                                             # site residue (47); union of legitimate
                                             # catalytic annotations only
        'source': 'UniProt Motif "Nudix box" (48-69) + Binding site residues (47, 49, 63, 67); '
                   'Suiter et al. 2020 PNAS',
        'reason': 'UniProt Region "Interaction with PCNA" (76-164) is a protein-protein '
                   'interaction site unrelated to catalytic function and was excluded; '
                   'it would otherwise have diluted the enrichment test by spanning 54% '
                   'of the protein.'
    },
    'PRKN': {
        'positions': [431, 238, 241, 253, 257, 260, 263, 289, 293, 332, 337, 352, 360,
                       365, 368, 373, 377, 418, 421, 436, 441, 446, 449, 457, 461],
        'source': 'UniProt Active site (431) + discrete Binding site residues, '
                   'Clausen et al. 2024 Nat Commun',
        'reason': 'UniProt Region annotations (77-237 "PINK1-dependent localization", '
                   '204-238/257-293 "SYT11 binding", 234-465 "TRIAD supradomain", '
                   '378-410 "REP") span up to 232 residues and are unrelated to direct '
                   'catalytic/ligase function; excluded to avoid diluting the enrichment '
                   'test (automated retrieval otherwise returned 389 of 465 positions, '
                   '84% of the protein).'
    },
    'VKORC1': {
        'positions': [43, 51, 55, 132, 135],
        'source': 'Catalytic C132XXC135 motif, C43/C51 cysteines involved in redox electron transfer, and '
                   'functionally constrained Phe55 hypothesized to participate in vitamin K binding;',
        'reason': 'UniProt automated retrieval returned only 3 sparse Binding site '
                   'positions (80, 135, 139), insufficient to represent the '
                   'well-characterized four-cysteine catalytic cycle plus substrate-binding '
                   'residue described in the primary literature.'
    }, 
    # Add overrides for other proteins ONLY if their UniProt auto-query proves sparse
    # OR contains annotations spanning unrelated biological functions 
}


AUTO_SITE_FEATURE_TYPES = ['Active site', 'Binding site', 'Site', 'Region', 'Motif']
# 'Domain' deliberately excluded - marks whole folded domains, not functional residues

MIN_AUTO_SITES_THRESHOLD = 5  # if UniProt returns fewer positions than this,
                                # flag for manual override consideration

AA3_TO_1 = {
    "Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E",
    "Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F",
    "Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V"
}

# ============================================================
# STEP 1: MaveDB download 
# ============================================================

def download_mavedb_scores(protein_name, urn):
    output_file = PROJECT_DIR / 'mavedb_data' / f'{protein_name}_scores.csv'
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    if output_file.exists():
        return output_file
    
    url = f"https://api.mavedb.org/api/v1/score-sets/{urn}/scores"
    response = requests.get(url, timeout=30)
    response.raise_for_status()
    with open(output_file, 'wb') as f:
        f.write(response.content)
    return output_file

# ============================================================
# STEP 2: AlphaMissense extraction 
# ============================================================

def extract_alphamissense(protein_name, uniprot_id, full_am_path=None):
    output_file = PROJECT_DIR / 'alphamissense_data' / f'AlphaMissense_{protein_name}.tsv'
    
    if output_file.exists():
        return pd.read_csv(output_file, sep='\t')
    
    if full_am_path is None:
        full_am_path = PROJECT_DIR / 'alphamissense_data' / 'AlphaMissense_aa_substitutions.tsv.gz'
    
    import gzip
    matched_rows, header = [], None
    with gzip.open(full_am_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            if header is None:
                header = line.strip().split('\t')
                continue
            if line.startswith(uniprot_id + '\t'):
                matched_rows.append(line.strip().split('\t'))
    
    df = pd.DataFrame(matched_rows, columns=header)
    df.to_csv(output_file, sep='\t', index=False)
    return df

# ============================================================
# STEP 3: Variant notation conversion 
# ============================================================

def hgvs3_to_1letter(hgvs_pro):
    if pd.isna(hgvs_pro):
        return None
    s = str(hgvs_pro).strip()
    if s.startswith("p."):
        s = s[2:]
    
    m = re.match(r'^([A-Z][a-z]{2})(\d+)=$', s)
    if m:
        wt3, pos = m.groups()
        return f"{AA3_TO_1.get(wt3,'?')}{pos}{AA3_TO_1.get(wt3,'?')}"
    m = re.match(r'^([A-Z][a-z]{2})(\d+)(Ter|\*)$', s)
    if m:
        wt3, pos, _ = m.groups()
        return f"{AA3_TO_1.get(wt3,'?')}{pos}*"
    m = re.match(r'^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$', s)
    if m:
        wt3, pos, mt3 = m.groups()
        return f"{AA3_TO_1.get(wt3,'?')}{pos}{AA3_TO_1.get(mt3,'?')}"
    return None

def a2a_transform(pathogenicity_score, a_min=0.05, a_max=0.95, k=6):
    p = np.asarray(pathogenicity_score, dtype=float)
    return a_min + (a_max - a_min) * (1.0 / (1.0 + np.exp(k * (p - 0.5))))

# ============================================================
# STEP 4: UniProt functional site query (with documented
#         fallback to manual override when sparse)
# ============================================================

def get_uniprot_sites(protein_name, uniprot_id):
    """
    Returns (final_positions, site_source, raw_auto_positions).
    raw_auto_positions is ALWAYS the unmodified automated UniProt result,
    even when a manual override is applied - needed to run the pure-UniProt
    sensitivity check regardless of curation decisions.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    response = requests.get(url)
    response.raise_for_status()
    data = response.json()
    
    raw_positions = set()
    print(f"  [{protein_name}] UniProt feature breakdown:")
    for feature in data.get('features', []):
        if feature.get('type') in AUTO_SITE_FEATURE_TYPES:
            loc = feature.get('location', {})
            start = loc.get('start', {}).get('value')
            end = loc.get('end', {}).get('value', start)
            if start is not None:
                span = end - start + 1
                desc = feature.get('description', '')
                print(f"    {feature.get('type')}: {start}-{end} "
                      f"({span} residues) — {desc}")
                if span > 20:
                    print(f"⚠️ SPAN >20 residues - likely too coarse to be a "
                          f"specific catalytic/binding site.")
                raw_positions.update(range(start, end + 1))
    
    print(f"  [{protein_name}] UniProt auto-detected {len(raw_positions)} site positions "
          f"(used for pure-UniProt sensitivity check regardless of override status)")

    if protein_name in MANUAL_SITE_OVERRIDES:
        override = MANUAL_SITE_OVERRIDES[protein_name]
        print(f"  [{protein_name}] Using documented manual override: {override['source']}")
        return set(override['positions']), 'manual_override', raw_positions
    
    if len(raw_positions) < MIN_AUTO_SITES_THRESHOLD:
        print(f"  [{protein_name}] WARNING: UniProt returned only {len(raw_positions)} sites, "
              f"no manual override defined. Enrichment test will likely lack power.")
        return raw_positions, 'uniprot_sparse', raw_positions
    
    return raw_positions, 'uniprot_auto', raw_positions

import requests
import pandas as pd

def export_uniprot_features(protein_name, uniprot_id):
    """
    Export complete UniProt feature retrieval for Supplementary Table S2.

    Returns
    -------
    list of dict
        One row per UniProt annotation.
    """

    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    response = requests.get(url)
    response.raise_for_status()
    data = response.json()

    rows = []

    manual_override = protein_name in MANUAL_SITE_OVERRIDES

    for feature in data.get("features", []):

        feature_type = feature.get("type")

        # Include Domain so reviewers can see why it was excluded
        if feature_type not in AUTO_SITE_FEATURE_TYPES and feature_type != "Domain":
            continue

        loc = feature.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value", start)

        if start is None:
            continue

        length = end - start + 1
        description = feature.get("description", "")

        decision = "Retained"
        reason = ""

        # ---------- DOMAIN ----------
        if feature_type == "Domain":
            decision = "Excluded"
            reason = "Broad structural domain"

        # ---------- PTEN ----------
        elif protein_name == "PTEN":
            if feature_type in ("Region", "Motif"):
                decision = "Excluded"
                reason = "Broad or non-catalytic annotation"
            else:
                decision = "Superseded"
                reason = "Manual site set replaced automated annotation"
        # ---------- NUDT15 ----------
        elif protein_name == "NUDT15":
            if feature_type == "Region" and "PCNA" in description:
                decision = "Excluded"
                reason = "Non-catalytic protein-protein interaction region"

        # ---------- PRKN ----------
        elif protein_name == "PRKN":
            if feature_type == "Region":
                decision = "Excluded"
                reason = "Broad functional region"

        # ---------- VKORC1 ----------
        elif protein_name == "VKORC1":
            decision = "Superseded"
            reason = "Manual site set replaced automated annotation"

        rows.append({
        "Protein": protein_name,
        "UniProt_ID": uniprot_id,
        "Feature_type": feature_type,
        "Start": start,
        "End": end,
        "Length": length,
        "Description": description,
        "Decision": decision,
        "Reason": reason,
        "Manual_override": "Yes" if manual_override else "No"
    })

    return rows

# ============================================================
# STEP 5: The complete per-protein pipeline (same function, 
#         called once for each protein in PROTEIN_CONFIGS)
# ============================================================

def run_threshold_sensitivity(valid, protein_name, functional_positions, score_column,
                                 pred_thresholds=(0.2, 0.3, 0.4), exp_thresholds=(0.7, 0.8, 0.9)):
    """
    Run enrichment test across a grid of failure-case threshold definitions.
    confirm this enrichment result is not an artifact of a
    single arbitrary threshold choice. Also computes the threshold-free
    continuous alternative (distance-to-site vs. residual correlation).
    """
    results = []
    for pred_t in pred_thresholds:
        for exp_t in exp_thresholds:
            is_failure = (valid['a2a_prediction'] < pred_t) & (valid[score_column] > exp_t)
            contingency = pd.crosstab(is_failure, valid['near_functional_site'])
            if contingency.shape == (2, 2):
                odds_ratio, p_val = stats.fisher_exact(contingency)
            else:
                odds_ratio, p_val = np.nan, np.nan
            results.append({
                'protein': protein_name, 'pred_threshold': pred_t, 'exp_threshold': exp_t,
                'n_failures': is_failure.sum(), 'odds_ratio': odds_ratio, 'p_value': p_val
            })
    
    sensitivity_df = pd.DataFrame(results)
    
    valid_ors = sensitivity_df['odds_ratio'].dropna()
    all_same_direction = (valid_ors > 1).all() or (valid_ors < 1).all()
    all_significant = (sensitivity_df['p_value'].dropna() < 0.05).all()
    
    print(f"  [{protein_name}] Threshold sensitivity: OR range "
          f"[{valid_ors.min():.2f}, {valid_ors.max():.2f}], "
          f"consistent direction: {all_same_direction}, "
          f"all significant: {all_significant}")
    
    def min_distance(pos):
        if pd.isna(pos) or not functional_positions:
            return np.nan
        return min(abs(pos - fp) for fp in functional_positions)
    
    valid = valid.copy()
    valid['dist_to_site'] = valid['position'].apply(min_distance)
    valid['residual'] = valid[score_column] - valid['a2a_prediction']
    
    if valid['dist_to_site'].notna().sum() > 10:
        cont_r, cont_p = stats.spearmanr(valid['dist_to_site'], valid['residual'])
        print(f"  [{protein_name}] Continuous check: distance-to-site vs residual "
              f"r={cont_r:.3f}, p={cont_p:.2e}")
    else:
        cont_r, cont_p = np.nan, np.nan
    
    return sensitivity_df, {
        'robust_direction': all_same_direction,
        'robust_significant': all_significant,
        'or_min': valid_ors.min() if len(valid_ors) else np.nan,
        'or_max': valid_ors.max() if len(valid_ors) else np.nan,
        'continuous_r': cont_r,
        'continuous_p': cont_p
    }

def odds_ratio_ci(contingency_table, alpha=0.05):
    """
    95% CI for an odds ratio using the Woolf/log-OR normal approximation.
    Standard method since scipy's fisher_exact doesn't provide a CI directly.
    Applies a 0.5 continuity correction if any cell is zero (Haldane-Anscombe),
    needed for PRKN's structural zinc-ligand subset where one cell is exactly 0.
    """
    table = contingency_table.values.astype(float)
    a, b = table[0, 0], table[0, 1]
    c, d = table[1, 0], table[1, 1]

    if 0 in (a, b, c, d):
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5

    log_or = np.log((a * d) / (b * c))
    se_log_or = np.sqrt(1/a + 1/b + 1/c + 1/d)
    z = stats.norm.ppf(1 - alpha / 2)

    ci_low = np.exp(log_or - z * se_log_or)
    ci_high = np.exp(log_or + z * se_log_or)
    return ci_low, ci_high

def run_pipeline_for_protein(protein_name, config):
    print(f"\n{'='*60}\nRunning pipeline for {protein_name}\n{'='*60}")
    
    # --- Load MaveDB ---
    mave_file = download_mavedb_scores(protein_name, config['mavedb_urn'])
    mavedb_df = pd.read_csv(mave_file)
    mavedb_df = mavedb_df[mavedb_df[config['hgvs_column']] != 'p.=']
    
    if config['score_column'] not in mavedb_df.columns:
        print(f"  ERROR: score_column '{config['score_column']}' not found.")
        print(f"  Available columns: {list(mavedb_df.columns)}")
        print(f"  --> Update PROTEIN_CONFIGS['{protein_name}']['score_column'] and rerun.")
        return None
    
    # --- Load AlphaMissense ---
    am_df = extract_alphamissense(protein_name, config['uniprot_id'])
    
    # --- Convert & merge ---
    mavedb_df['variant_1letter'] = mavedb_df[config['hgvs_column']].apply(hgvs3_to_1letter)
    merged = mavedb_df.merge(
        am_df[['protein_variant', 'am_pathogenicity', 'am_class']],
        left_on='variant_1letter', right_on='protein_variant', how='inner'
    )
    coverage = 100 * len(merged) / len(mavedb_df)
    print(f"  Matched: {len(merged)}/{len(mavedb_df)} ({coverage:.1f}%)")

    # Force numeric dtype - fresh AlphaMissense extractions (built from raw split text lines) 
    merged['am_pathogenicity'] = pd.to_numeric(merged['am_pathogenicity'], errors='coerce')
    n_bad = merged['am_pathogenicity'].isna().sum()
    if n_bad > 0:
        print(f"  WARNING: {n_bad} rows had non-numeric am_pathogenicity values, dropped")
        merged = merged.dropna(subset=['am_pathogenicity'])
    
    # --- A2A transform & correlation ---
    merged['a2a_prediction'] = a2a_transform(merged['am_pathogenicity'])
    valid = merged.dropna(subset=[config['score_column'], 'a2a_prediction'])
    
    spearman_r, spearman_p = stats.spearmanr(valid['a2a_prediction'], valid[config['score_column']])
    print(f"  Spearman r = {spearman_r:.3f} (p={spearman_p:.2e}), n={len(valid)}")
    
    # --- Baseline check: does the A2A transformation add anything over raw ---
    # --- AlphaMissense pathogenicity, simply inverted? ---
    from sklearn.metrics import mean_absolute_error
    raw_inverted = 1.0 - valid['am_pathogenicity']

    pearson_a2a, pearson_a2a_p = stats.pearsonr(valid['a2a_prediction'], valid[config['score_column']])
    pearson_raw, pearson_raw_p = stats.pearsonr(raw_inverted, valid[config['score_column']])

    mae_a2a = mean_absolute_error(valid[config['score_column']], valid['a2a_prediction'])
    mae_raw = mean_absolute_error(valid[config['score_column']], raw_inverted)

    print(f"  Calibration check (Pearson r / MAE vs experimental score):")
    print(f"    A2A:          r={pearson_a2a:.3f} (p={pearson_a2a_p:.2e}), MAE={mae_a2a:.3f}")
    print(f"    Raw (1-P):    r={pearson_raw:.3f} (p={pearson_raw_p:.2e}), MAE={mae_raw:.3f}")
    print(f"    Delta: Pearson r {pearson_a2a - pearson_raw:+.3f}, MAE {mae_a2a - mae_raw:+.3f} "
          f"(negative MAE delta = A2A more accurate)")

    # --- CRITICAL FOLLOW-UP CHECK ---
    valid['raw_inverted'] = raw_inverted
    valid['is_failure_case_raw'] = (
        (valid['raw_inverted'] < 0.4) & (valid[config['score_column']] > 0.8)
    )
    n_diff_classification = (valid['is_failure_case_raw'] != 
                              ((valid['a2a_prediction'] < 0.4) & (valid[config['score_column']] > 0.8))).sum()
    print(f"    Variants classified differently as failure-case under raw vs A2A: "
          f"{n_diff_classification}/{len(valid)} ({100*n_diff_classification/len(valid):.1f}%)")
    
    # --- UniProt site enrichment ---
    def extract_position(v):
        m = re.match(r'^[A-Z](\d+)[A-Z]$', str(v))
        return int(m.group(1)) if m else None
    valid['position'] = valid['variant_1letter'].apply(extract_position)
    
    functional_positions, site_source, raw_uniprot_positions = get_uniprot_sites(protein_name, config['uniprot_id'])
    valid['near_functional_site'] = valid['position'].apply(
        lambda p: p in functional_positions if pd.notna(p) else False
    )
    
    valid['is_failure_case'] = (
        (valid['a2a_prediction'] < 0.4) & (valid[config['score_column']] > 0.8)
    )
    
    contingency = pd.crosstab(valid['is_failure_case'], valid['near_functional_site'])
    if contingency.shape == (2, 2) and valid['near_functional_site'].sum() >= 5:
        odds_ratio, p_value = stats.fisher_exact(contingency)
        ci_low, ci_high = odds_ratio_ci(contingency)
        print(f"  Enrichment (A2A-based): OR={odds_ratio:.2f} [95% CI: {ci_low:.2f}-{ci_high:.2f}], "
              f"p={p_value:.2e} (site source: {site_source})")
    else:
        odds_ratio, p_value, ci_low, ci_high = np.nan, np.nan, np.nan, np.nan
        print(f"  Enrichment test skipped - insufficient sites (source: {site_source})")

    # --- Sensitivity check: pure UniProt annotations, no manual curation ---
    # --- Run ONLY when a manual override was applied, to test whether the ---
    # --- qualitative conclusion depends on curation, or holds even with ---
    # --- unmodified (potentially imprecise or diluted) automated annotations. ---
    pure_or, pure_p, pure_n_sites = np.nan, np.nan, len(raw_uniprot_positions)
    if site_source == 'manual_override':
        valid['near_raw_uniprot_site'] = valid['position'].apply(
            lambda p: p in raw_uniprot_positions if pd.notna(p) else False
        )
        pure_contingency = pd.crosstab(valid['is_failure_case'], valid['near_raw_uniprot_site'])
        if pure_contingency.shape == (2, 2) and valid['near_raw_uniprot_site'].sum() >= 5:
            pure_or, pure_p = stats.fisher_exact(pure_contingency)
            print(f"  SENSITIVITY (pure UniProt, {pure_n_sites} sites, no curation): "
                  f"OR={pure_or:.2f}, p={pure_p:.2e}")
        else:
            print(f"  SENSITIVITY (pure UniProt, {pure_n_sites} sites, no curation): "
                  f"insufficient sites to test (n_near={valid['near_raw_uniprot_site'].sum()})")

    # --- Robustness check: does enrichment hold using RAW (1-P) classification ---
    # --- instead of A2A? This directly tests whether A2A's specific transformation ---
    # --- is necessary for the enrichment finding, or whether raw pathogenicity alone suffices. ---
    contingency_raw = pd.crosstab(valid['is_failure_case_raw'], valid['near_functional_site'])
    if contingency_raw.shape == (2, 2) and valid['near_functional_site'].sum() >= 5:
        odds_ratio_raw, p_value_raw = stats.fisher_exact(contingency_raw)
        print(f"  Enrichment (raw 1-P based): OR={odds_ratio_raw:.2f}, p={p_value_raw:.2e}")
    else:
        odds_ratio_raw, p_value_raw = np.nan, np.nan
        print(f"  Enrichment (raw 1-P based): test skipped - insufficient sites")
    
    # --- Threshold sensitivity + continuous robustness check (standard for every protein) ---
    sensitivity_df, robustness = run_threshold_sensitivity(
        valid, protein_name, functional_positions, config['score_column']
    )
    sensitivity_df.to_csv(PROJECT_DIR / f'{protein_name}_threshold_sensitivity.csv', index=False)
    
    # --- Save results ---
    out_file = PROJECT_DIR / f'{protein_name}_final_results.csv'
    valid.to_csv(out_file, index=False)
    
    return {
        'protein': protein_name, 'n_matched': len(valid), 'coverage_pct': coverage,
        'spearman_r': spearman_r, 'spearman_p': spearman_p,
        'pearson_a2a': pearson_a2a, 'pearson_raw': pearson_raw,
        'mae_a2a': mae_a2a, 'mae_raw': mae_raw,
        'enrichment_or': odds_ratio, 'enrichment_p': p_value,
        'enrichment_ci_low': ci_low, 'enrichment_ci_high': ci_high,
        'pure_uniprot_or': pure_or, 'pure_uniprot_p': pure_p, 'pure_uniprot_n_sites': pure_n_sites,
        'enrichment_or_raw': odds_ratio_raw, 'enrichment_p_raw': p_value_raw,
        'site_source': site_source,
        'robust_direction': robustness['robust_direction'],
        'robust_significant': robustness['robust_significant'],
        'or_range': f"[{robustness['or_min']:.2f}, {robustness['or_max']:.2f}]",
        'continuous_r': robustness['continuous_r'],
        'continuous_p': robustness['continuous_p'],
    }

# ============================================================
# MAIN — run for all configured proteins, one summary table
# ============================================================

if __name__ == "__main__":
    all_results = []
    supplementary_rows = []

    for protein_name, config in PROTEIN_CONFIGS.items():
        result = run_pipeline_for_protein(protein_name, config)
        if result:
            all_results.append(result)
    

        supplementary_rows.extend(
            export_uniprot_features(protein_name, config['uniprot_id'])
        )

    summary_df = pd.DataFrame(all_results)

    # ============================================================
    # Multiple-testing correction (Benjamini-Hochberg) across the
    # six main enrichment tests 
    # ============================================================
    from statsmodels.stats.multitest import multipletests

    valid_p_mask = summary_df['enrichment_p'].notna()
    if valid_p_mask.sum() > 0:
        rejected, p_adj, _, _ = multipletests(
            summary_df.loc[valid_p_mask, 'enrichment_p'], method='fdr_bh'
        )
        summary_df.loc[valid_p_mask, 'enrichment_p_fdr'] = p_adj
        summary_df.loc[valid_p_mask, 'enrichment_significant_fdr'] = rejected

        print(f"\n{'='*60}")
        print("BENJAMINI-HOCHBERG FDR CORRECTION (main enrichment tests)")
        print(f"{'='*60}")
        for _, row in summary_df.loc[valid_p_mask].iterrows():
            print(f"  {row['protein']}: raw p={row['enrichment_p']:.2e} -> "
                  f"FDR-adjusted p={row['enrichment_p_fdr']:.2e} "
                  f"(significant: {row['enrichment_significant_fdr']})")

    summary_df.to_csv(PROJECT_DIR / 'all_proteins_summary.csv', index=False)
    
    print(f"\n{'='*60}\nSUMMARY ACROSS ALL PROTEINS\n{'='*60}")
    print(summary_df.to_string(index=False))
    #--------------------
    
    supp_df = pd.DataFrame(supplementary_rows)

    # Remove duplicated UniProt annotations
    supp_df = supp_df.drop_duplicates()

    supp_df = supp_df.sort_values(
        by=["Protein", "Feature_type", "Start"]
    )

    supp_df.to_csv("pathogenicity_abundance_project/Supplementary_Table_S1.csv", index=False)

    print("Saved Supplementary_Table_S1.csv")

