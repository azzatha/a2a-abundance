"""
Data Download for Pathogenicity-Abundance Project — FINAL VERSION
==================================================================
Downloads/prepares:
  1. MaveDB variant abundance scores for all six proteins
  2. AlphaMissense predictions (via full-file extraction)
"""

import gzip
import requests
import pandas as pd
from pathlib import Path

# ============================================================
# PATHS — previously missing, caused a NameError on import
# ============================================================
PROJECT_DIR = Path('pathogenicity_abundance_project')
MAVEDB_DIR = PROJECT_DIR / 'mavedb_data'
ALPHAMISSENSE_DIR = PROJECT_DIR / 'alphamissense_data'
MAVEDB_DIR.mkdir(parents=True, exist_ok=True)
ALPHAMISSENSE_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# PART 1: MaveDB Downloads — all six proteins, all URNs verified
# ============================================================
DATASETS = {
    'PTEN': {
        'urn': 'urn:mavedb:00000013-a-1',
        'uniprot': 'P60484',
        'verified': True,
        'note': 'Matreyek et al. 2018, Nat Genet.'
    },
    'TPMT': {
        'urn': 'urn:mavedb:00000013-b-1',
        'uniprot': 'P51580',
        'verified': True,
        'note': 'Matreyek et al. 2018, Nat Genet.'
    },
    'NUDT15': {
        'urn': 'urn:mavedb:00000055-a-1',
        'uniprot': 'Q9NV35',
        'verified': True,
        'note': 'Suiter et al. 2020, PNAS.'
    },
    'ASPA': {
        'urn': 'urn:mavedb:00000657-a-1',
        'uniprot': 'P45381',
        'verified': True,
        'note': 'Gronbaek-Thygesen et al. 2024, Nat Commun.'
    },
    'PRKN': {
        'urn': 'urn:mavedb:00000114-a-1',
        'uniprot': 'O60260',
        'verified': True,
        'note': 'Clausen et al. 2024, Nat Commun. Original abundance dataset '
                'matching Livesey & Marsh 2025 (NOT the newer 2026 '
                'depolarisation-specific set, urn:mavedb:00001281-a-1).'
    },
    'VKORC1': {
        'urn': 'urn:mavedb:00000078-b-1',
        'uniprot': 'Q9BQB6',
        'verified': True,
        'note': 'Chiasson et al. 2020, eLife.'
    },
}

MAVEDB_API_BASE = "https://api.mavedb.org/api/v1"

def download_mavedb_scores(protein_name, urn, output_dir=MAVEDB_DIR):
    output_file = output_dir / f"{protein_name}_scores.csv"
    if output_file.exists():
        print(f"  [{protein_name}] Already downloaded -> {output_file}")
        return output_file

    url = f"{MAVEDB_API_BASE}/score-sets/{urn}/scores"
    print(f"  [{protein_name}] Fetching {url} ...")
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        with open(output_file, 'wb') as f:
            f.write(response.content)
        print(f"  [{protein_name}] Saved -> {output_file}")
        return output_file
    except Exception as e:
        print(f"  [{protein_name}] API FAILED ({e})")
        print(f"  [{protein_name}] Manual fallback: visit "
              f"https://www.mavedb.org/score-sets/{urn} and click 'Download scores',")
        print(f"  [{protein_name}] then save the file to: {output_file}")
        return None

def download_all_mavedb():
    print("\n" + "="*60 + "\nPART 1: Downloading MaveDB scores\n" + "="*60)
    results = {}
    for protein, info in DATASETS.items():
        if not info.get('verified', False):
            print(f"  [{protein}] SKIPPED - URN not yet verified.")
            continue
        results[protein] = download_mavedb_scores(protein, info['urn'])
    return results

# ============================================================
# PART 2: AlphaMissense Downloads
# ============================================================
ALPHAMISSENSE_URL = "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz"
ALPHAMISSENSE_FULL_PATH = ALPHAMISSENSE_DIR / "AlphaMissense_aa_substitutions.tsv.gz"

def download_full_alphamissense():
    if ALPHAMISSENSE_FULL_PATH.exists():
        size_mb = ALPHAMISSENSE_FULL_PATH.stat().st_size / 1e6
        print(f"Already downloaded ({size_mb:.0f} MB) -> {ALPHAMISSENSE_FULL_PATH}")
        return ALPHAMISSENSE_FULL_PATH

    print(f"Downloading AlphaMissense full dataset from:\n  {ALPHAMISSENSE_URL}")
    try:
        response = requests.get(ALPHAMISSENSE_URL, stream=True, timeout=60)
        response.raise_for_status()
        total = int(response.headers.get('content-length', 0))
        downloaded = 0
        with open(ALPHAMISSENSE_FULL_PATH, 'wb') as f:
            for chunk in response.iter_content(chunk_size=1024*1024):
                f.write(chunk)
                downloaded += len(chunk)
                if total > 0:
                    print(f"\r  Progress: {100*downloaded/total:.1f}%", end='')
        print(f"\nDone -> {ALPHAMISSENSE_FULL_PATH}")
        return ALPHAMISSENSE_FULL_PATH
    except Exception as e:
        print(f"\nDOWNLOAD FAILED: {e}")
        print("Manual fallback: github.com/google-deepmind/alphamissense")
        return None

def extract_protein_from_alphamissense(uniprot_id, protein_name, full_file_path=ALPHAMISSENSE_FULL_PATH):
    output_file = ALPHAMISSENSE_DIR / f"AlphaMissense_{protein_name}.tsv"
    if output_file.exists():
        print(f"  [{protein_name}] Already extracted -> {output_file}")
        return pd.read_csv(output_file, sep='\t')

    if not full_file_path or not Path(full_file_path).exists():
        print(f"  [{protein_name}] SKIPPED - full AlphaMissense file not found.")
        return None

    print(f"  [{protein_name}] Scanning full file for UniProt ID {uniprot_id} ...")
    matched_rows, header = [], None
    with gzip.open(full_file_path, 'rt') as f:
        for i, line in enumerate(f):
            if line.startswith('#'):
                continue
            if header is None:
                header = line.strip().split('\t')
                continue
            if line.startswith(uniprot_id + '\t'):
                matched_rows.append(line.strip().split('\t'))
            if i % 10_000_000 == 0 and i > 0:
                print(f"    ...scanned {i:,} lines, {len(matched_rows)} matches")

    if not matched_rows:
        print(f"  [{protein_name}] WARNING: no rows found for {uniprot_id}.")
        return None

    df = pd.DataFrame(matched_rows, columns=header)
    df.to_csv(output_file, sep='\t', index=False)
    print(f"  [{protein_name}] Extracted {len(df)} rows -> {output_file}")
    return df

def download_all_alphamissense():
    print("\n" + "="*60 + "\nPART 2: AlphaMissense predictions\n" + "="*60)
    full_file = download_full_alphamissense()
    results = {}
    for protein, info in DATASETS.items():
        if not info.get('verified', False):
            continue
        results[protein] = extract_protein_from_alphamissense(info['uniprot'], protein, full_file)
    return results

if __name__ == "__main__":
    print("Starting data download pipeline...\n")
    mavedb_results = download_all_mavedb()
    am_results = download_all_alphamissense()

    print("\n" + "="*60 + "\nDOWNLOAD SUMMARY\n" + "="*60)
    for protein in DATASETS:
        mave_ok = mavedb_results.get(protein) is not None
        am_ok = am_results.get(protein) is not None
        print(f"  {protein}: MaveDB={'OK' if mave_ok else 'FAILED'}, "
              f"AlphaMissense={'OK' if am_ok else 'FAILED'}")
