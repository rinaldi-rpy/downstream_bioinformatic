#!/usr/bin/env bash
# Select stages explicitly instead of commenting out executable code.
set -Eeuo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
if [[ ${1:-} == --help ]]; then
    echo 'Usage: STEPS=qc,ld,pca,admixture,fst,bootstrap,treemix bash run_pipeline.sh INPUT.vcf.gz METADATA.csv'
    echo 'Choose any subset with STEPS. Default: qc,ld,pca,admixture,fst. See README.md.'
    exit 0
fi
source "$SCRIPT_DIR/shell/common.sh"
VCF=${1:?Input VCF required}; META=${2:?Metadata required}
VCF=$(realpath "$VCF"); META=$(realpath "$META")
STEPS=${STEPS:-qc,ld,pca,admixture,fst}
RUN_DIR=${RUN_DIR:-$BASE_DIR/runs/analysis_$(date +%Y%m%d_%H%M%S)_$$}
# INPUT_RUN supplies pre-existing stage products for partial reruns.
INPUT_RUN=${INPUT_RUN:-$RUN_DIR}
[[ ! -e "$RUN_DIR" ]] || { echo 'RUN_DIR must be a fresh directory.' >&2; exit 1; }
require_file "$VCF"; require_file "$META"; require_file "$AUTOSOME_MAP"
"$PYTHON" "$SCRIPT_DIR/python/metadata.py" summary --meta "$META" --group "$GROUP_COLUMN"
IFS=, read -r -a stages <<< "$STEPS"
for stage in "${stages[@]}"; do
    case "$stage" in qc|ld|pca|admixture|fst|bootstrap|treemix) ;; *) echo "Unknown stage: $stage" >&2; exit 1;; esac
done
mkdir -p "$RUN_DIR"; RUN_DIR=$(realpath "$RUN_DIR"); INPUT_RUN=$(realpath -m "$INPUT_RUN")
exec > >(tee "$RUN_DIR/pipeline.log") 2>&1
printf '%s\n' "VCF=$VCF" "META=$META" "INPUT_RUN=$INPUT_RUN" "STEPS=$STEPS" > "$RUN_DIR/inputs.txt"
cp "$META" "$RUN_DIR/metadata.${META##*.}"; cp "$AUTOSOME_MAP" "$RUN_DIR/autosomes.tsv"
# Record resolved parameters and source hashes without dumping unrelated environment values.
"$PYTHON" "$SCRIPT_DIR/python/provenance.py" "$SCRIPT_DIR" "$RUN_DIR"
QC_VCF="$INPUT_RUN/01_qc/03.filtered.recode.vcf.gz"
LD_PREFIX="$INPUT_RUN/02_ld/merged.pruned"
AUTO_VCF="$INPUT_RUN/02_ld/autosomal.vcf.gz"
for stage in "${stages[@]}"; do
    echo "Running stage: $stage"
    case "$stage" in
        qc) source "$SCRIPT_DIR/shell/qc.sh"; QC_VCF="$RUN_DIR/01_qc/03.filtered.recode.vcf.gz" ;;
        ld) source "$SCRIPT_DIR/shell/ld.sh"; LD_PREFIX="$RUN_DIR/02_ld/merged.pruned"; AUTO_VCF="$RUN_DIR/02_ld/autosomal.vcf.gz" ;;
        pca) source "$SCRIPT_DIR/shell/pca.sh" ;;
        admixture) source "$SCRIPT_DIR/shell/admixture.sh" ;;
        fst) source "$SCRIPT_DIR/shell/fst.sh" ;;
        bootstrap) BASE_DIR="$BASE_DIR" VCF_PATH="$AUTO_VCF" META_PATH="$META" OUT_DIR="$RUN_DIR/05_bootstrap" bash "$SCRIPT_DIR/run_fst_bootstrap.sh" ;;
        treemix) BASE_DIR="$BASE_DIR" BFILE="$LD_PREFIX" META_PATH="$META" OUT_DIR="$RUN_DIR/06_treemix" bash "$SCRIPT_DIR/run_treemix.sh" ;;
    esac
done
printf 'Completed requested stages. Results: %s\n' "$RUN_DIR"
