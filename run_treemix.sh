#!/usr/bin/env bash
# Modular TreeMix pipeline for the buffalo dataset.
# Required: Bash, Python 3, PLINK 1.9, TreeMix.
# Optional plots: Rscript and RColorBrewer. The local TreeMix helper is versioned with this workflow.
# Run: bash run_treemix.sh
# Override settings: ROOT=BGD-NE MAX_M=3 REPS=10 BLOCK_SIZE=500 bash run_treemix.sh
# Repeated random seeds assess optimization stability; these are not bootstraps.
set -Eeuo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/shell/common.sh"
trap 'echo "ERROR at line $LINENO. Check the pipeline and per-run logs." >&2' ERR

BASE_DIR="${BASE_DIR:-$SCRIPT_DIR}"
BFILE="${BFILE:-$BASE_DIR/work/02_ld/merged.pruned}"
META_PATH="${META_PATH:-$BASE_DIR/metadata.csv}"
# Use a fresh directory by default to keep results from different settings separate.
OUT_DIR="${OUT_DIR:-$BASE_DIR/work/05_treemix/run_$(date +%Y%m%d_%H%M%S)_$$}"
MAX_M="${MAX_M:-8}"
REPS="${REPS:-10}"
BLOCK_SIZE="${BLOCK_SIZE:-500}"
ROOT="${ROOT:-}"
PLINK="${PLINK:-plink}"
TREEMIX="${TREEMIX:-treemix}"
PYTHON="${PYTHON:-python3}"
PLOT_FUNCS="${PLOT_FUNCS:-}"
unset RUN_DIR
export AUTOSOME_MAP GROUP_COLUMN BASE_DIR BFILE META_PATH OUT_DIR MAX_M REPS BLOCK_SIZE ROOT PLOT_FUNCS

for executable in "$PLINK" "$TREEMIX" "$PYTHON"; do
    command -v "$executable" >/dev/null || { echo "Missing executable: $executable" >&2; exit 1; }
done
for suffix in bed bim fam; do
    [[ -s "$BFILE.$suffix" ]] || { echo "Missing input: $BFILE.$suffix" >&2; exit 1; }
done
[[ -s "$META_PATH" ]] || { echo "Missing metadata: $META_PATH" >&2; exit 1; }
[[ "$MAX_M" =~ ^[0-9]+$ && "$REPS" =~ ^[1-9][0-9]*$ && "$BLOCK_SIZE" =~ ^[1-9][0-9]*$ ]] || {
    echo "MAX_M must be nonnegative; REPS and BLOCK_SIZE must be positive integers." >&2; exit 1;
}
[[ ! -e "$OUT_DIR" ]] || { echo "Output directory already exists; choose a new OUT_DIR." >&2; exit 1; }
mkdir -p "$OUT_DIR" "$OUT_DIR/runs" "$OUT_DIR/best" "$OUT_DIR/plots"
exec > >(tee "$OUT_DIR/pipeline.log") 2>&1
plink_version=$("$PLINK" --version)
[[ "$plink_version" == *"PLINK v1.9"* ]] || { echo "This pipeline requires PLINK 1.9: $plink_version"; exit 1; }
printf '%s\n' "$plink_version"
printf 'BFILE=%s\nMETA_PATH=%s\nMAX_M=%s\nREPS=%s\nBLOCK_SIZE=%s\nROOT=%s\n' \
    "$BFILE" "$META_PATH" "$MAX_M" "$REPS" "$BLOCK_SIZE" "$ROOT" > "$OUT_DIR/settings.txt"

"$PYTHON" "$SCRIPT_DIR/python/provenance.py" "$SCRIPT_DIR" "$OUT_DIR"
echo "Step 1: Match metadata and select autosomal SNPs."
"$PYTHON" "${SCRIPT_DIR}/python/treemix_map.py"

echo "Step 2: Calculate population-stratified allele counts."
"$PLINK" --bfile "$BFILE" --chr-set 24 no-xy --allow-extra-chr --keep-allele-order \
    --extract "$OUT_DIR/autosomal_snps.txt" --nonfounders \
    --within "$OUT_DIR/populations.txt" --freq --out "$OUT_DIR/autosomal"

echo "Step 3: Convert allele counts to TreeMix input."
"$PYTHON" "${SCRIPT_DIR}/python/treemix_convert.py"

echo "Step 4: Run TreeMix models sequentially with independent random seeds."
printf "M\tReplicate\tSeed\n" > "$OUT_DIR/seeds.tsv"
root_args=()
[[ -z "$ROOT" ]] || root_args=(-root "$ROOT")
for ((m=0; m<=MAX_M; m++)); do
    for ((rep=1; rep<=REPS; rep++)); do
        seed=$((SEED + m * REPS + rep))
        printf "%s\t%s\t%s\n" "$m" "$rep" "$seed" >> "$OUT_DIR/seeds.tsv"
        prefix="$OUT_DIR/runs/treemix_m${m}_rep${rep}"
        echo "Running m=$m, replicate=$rep, seed=$seed"
        "$TREEMIX" -i "$OUT_DIR/treemix_input.gz" -m "$m" -k "$BLOCK_SIZE" \
            -seed "$seed" -global "${root_args[@]}" -o "$prefix" > "$prefix.log" 2>&1
        [[ -s "$prefix.llik" && -s "$prefix.treeout.gz" ]] || {
            echo "Incomplete output: $prefix; inspect its log."; exit 1;
        }
    done
done

echo "Step 5: Summarize likelihoods and select the best replicate per model."
"$PYTHON" "${SCRIPT_DIR}/python/treemix_summary.py"

echo "Step 6: Generate plots when plotting dependencies are available."
PLOT_FUNCS="${PLOT_FUNCS:-$SCRIPT_DIR/r/treemix_plot_functions.R}"
export PLOT_FUNCS
if command -v Rscript >/dev/null; then
    Rscript --vanilla "$SCRIPT_DIR/r/plot_treemix.R"
else
    echo "Plots skipped: Rscript unavailable. Numeric outputs are retained."
fi
echo "Pipeline completed. Output directory: $OUT_DIR"
echo "Review model_summary.csv, replicate stability, and residuals before choosing migration count."
echo "This run does not estimate bootstrap support or establish migration direction without rooting review."
