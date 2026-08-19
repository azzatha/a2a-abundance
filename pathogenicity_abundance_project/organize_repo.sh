#!/usr/bin/env bash
# Build the a2a-abundance repository structure by COPYING from the existing
# project folder. Nothing is moved or deleted, so your current layout is
# untouched and you can re-run this safely.
#
# Adjust the two paths below, then run:  bash organize_repo.sh

set -u

SRC="pathogenicity_abundance_project"   # existing project folder
CODE_SRC="."                            # where your .py scripts live
DEST="a2a-abundance"                    # new repo folder (created)

PROTEINS=(PTEN TPMT NUDT15 ASPA PRKN VKORC1)

echo "Building $DEST/ ..."
mkdir -p "$DEST"/{code,data/mavedb,data/alphamissense,results,figures}

copy () {  # copy $1 -> $2, report if missing
  if [ -e "$1" ]; then cp -r "$1" "$2" && echo "  ok    $1"
  else echo "  MISS  $1"; fi
}

echo
echo "[code]"
copy "$CODE_SRC/unified_pipeline-final.py" "$DEST/code/unified_pipeline.py"
copy "$CODE_SRC/prkn_site_isolation.py"    "$DEST/code/"
copy "$CODE_SRC/make_all_figures.py"       "$DEST/code/"
copy "$CODE_SRC/verify_prkn.py"            "$DEST/code/"

echo
echo "[results]"
copy "$SRC/all_proteins_summary.csv"       "$DEST/results/"
copy "$SRC/Supplementary_Table_S1.csv"     "$DEST/results/"
for p in "${PROTEINS[@]}"; do
  copy "$SRC/${p}_final_results.csv"        "$DEST/results/"
  copy "$SRC/${p}_threshold_sensitivity.csv" "$DEST/results/"
done

echo
echo "[input data - committed so results reproduce without network access]"
copy "$SRC/mavedb_data/."                  "$DEST/data/mavedb/"
# per-protein AlphaMissense extracts only; the ~1 GB source file is excluded
for p in "${PROTEINS[@]}"; do
  copy "$SRC/alphamissense_data/AlphaMissense_${p}.tsv" "$DEST/data/alphamissense/"
done

echo
echo "[figures]"
for f in figure1_scatter figure2_prkn_decomposition figureS1_threshold_heatmap; do
  copy "figures/${f}.pdf" "$DEST/figures/"
  copy "figures/${f}.png" "$DEST/figures/"
done

echo
echo "[repo files]"
copy "README.md"       "$DEST/"
copy "requirements.txt" "$DEST/"
copy ".gitignore"      "$DEST/"

echo
echo "Checking nothing oversized slipped in (GitHub limit is 100 MB/file):"
find "$DEST" -type f -size +50M -exec ls -lh {} \; | awk '{print "  LARGE  " $9 " (" $5 ")"}'
echo "  (no output above = all files fine)"

echo
echo "Total size:"; du -sh "$DEST"
echo "Done. Review $DEST/, then follow the git steps."
