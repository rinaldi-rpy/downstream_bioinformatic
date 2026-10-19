#!/usr/bin/env bash
# Resume TreeMix postprocessing only. Does not run PLINK or TreeMix.
# Compatible with Python 3.9+ (no assignment expressions).
# Usage: bash summarize_treemix.sh /path/to/existing/run_directory
set -Eeuo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/shell/common.sh"
trap 'echo "ERROR at line $LINENO. Original TreeMix outputs are unchanged." >&2' ERR
PYTHON="${PYTHON:-python3}"
PLOT_FUNCS="${PLOT_FUNCS:-}"
RUN_DIR="${1:-${RUN_DIR:-}}"
SEARCH_DIR="${SEARCH_DIR:-$SCRIPT_DIR/work/05_treemix}"
export RUN_DIR SEARCH_DIR
command -v "$PYTHON" >/dev/null || { echo "Python 3 is required."; exit 1; }
RUN_DIR=$("$PYTHON" "${SCRIPT_DIR}/python/treemix_find_run.py"
)
export RUN_DIR
# Read original settings; discover run counts only if settings were not saved.
settings=$("$PYTHON" "${SCRIPT_DIR}/python/treemix_read_settings.py"
)
read -r MAX_M REPS <<< "$settings"
OUT_DIR="$RUN_DIR/postprocess_$(date +%Y%m%d_%H%M%S)_$$"
export OUT_DIR MAX_M REPS PLOT_FUNCS
mkdir -p "$OUT_DIR/best" "$OUT_DIR/plots"
exec > >(tee "$OUT_DIR/postprocess.log") 2>&1
echo "Existing run: $RUN_DIR"
echo "Postprocessing output: $OUT_DIR"
echo "Models: m=0..$MAX_M; replicates per model: $REPS"
"$PYTHON" "${SCRIPT_DIR}/python/treemix_copy_population_order.py"
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
echo "Postprocessing completed. Output directory: $OUT_DIR"
echo "Review model_summary.csv, replicate stability, and residuals before choosing migration count."
echo "This run does not estimate bootstrap support or establish migration direction without rooting review."
