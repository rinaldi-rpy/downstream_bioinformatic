#!/usr/bin/env bash
# Streaming descriptive QC; paths are explicit and INFO/DP is distinct from FORMAT/DP.
set -Eeuo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/shell/common.sh"
[[ $# == 4 ]] || { echo 'Usage: bash summarize_qc.sh RAW.vcf.gz QC.vcf.gz PRUNED.bim OUTPUT.csv' >&2; exit 1; }
"$PYTHON" "$SCRIPT_DIR/python/qc_summary.py" --raw "$1" --qc "$2" --bim "$3" --out "$4"
