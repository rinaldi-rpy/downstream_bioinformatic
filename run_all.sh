#!/usr/bin/env bash
# One command for the complete revised buffalo workflow.
# This launcher changes stage selection only, not scientific parameters.
set -Eeuo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
usage() {
    echo 'Usage: bash run_all.sh INPUT.vcf.gz METADATA.csv'
    echo 'Stages: qc -> ld -> pca -> admixture -> fst -> bootstrap -> treemix'
    echo 'Configuration: config/parameters.sh; environment overrides are supported.'
    echo 'Results: a fresh runs/analysis_<timestamp>_<pid>/ directory, or RUN_DIR.'
}
if [[ ${1:-} == --help || ${1:-} == -h ]]; then
    usage
    exit 0
fi
if [[ $# -ne 2 ]]; then
    usage >&2
    exit 2
fi
source "$SCRIPT_DIR/shell/common.sh"
require_file "$1"
require_file "$2"
# Fail before expensive QC if a later mandatory tool/package is unavailable.
for executable in "$BCFTOOLS" "$VCFTOOLS" "$ADMIXTURE" "$TREEMIX" Rscript; do
    need "$executable"
done
require_plink19
"$PYTHON" -c 'import numpy; print("NumPy:", numpy.__version__)'
Rscript --vanilla -e 'for (p in c("ape", "RColorBrewer")) if (!requireNamespace(p, quietly=TRUE)) stop(paste("Missing R package:", p))'
export STEPS=qc,ld,pca,admixture,fst,bootstrap,treemix
exec bash "$SCRIPT_DIR/run_pipeline.sh" "$1" "$2"
