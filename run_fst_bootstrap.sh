#!/usr/bin/env bash
# Pairwise Weir-Cockerham FST, shared genomic block bootstrap, and NJ tree.
# Required: Bash, PLINK 1.9, Python >=3.9 with NumPy >=1.17, R with ape.
# Default input is the unpruned QC VCF, not the overlapping window FST tables.
# Usage: bash run_fst_bootstrap.sh
# Example: N_BOOT=2000 BLOCK_BP=1000000 bash run_fst_bootstrap.sh
# Optional BFILE=/path/to/unpruned/prefix bypasses VCF conversion.
set -Eeuo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/shell/common.sh"
trap 'echo "ERROR at line $LINENO. Check pipeline.log; original data are unchanged." >&2' ERR
BASE_DIR="${BASE_DIR:-$SCRIPT_DIR}"
VCF_PATH="${VCF_PATH:-$BASE_DIR/work/01_qc/03.filtered.recode.vcf.gz}"
META_PATH="${META_PATH:-$BASE_DIR/metadata.csv}"
FST_DIR="${FST_DIR:-$BASE_DIR/work/05_fst}"
OUT_DIR="${OUT_DIR:-$FST_DIR/bootstrap_nj_$(date +%Y%m%d_%H%M%S)_$$}"
BFILE="${BFILE:-}"
BLOCK_BP="${BLOCK_BP:-1000000}"
N_BOOT="${N_BOOT:-2000}"
SEED="${SEED:-20261015}"
CHUNK_SNPS="${CHUNK_SNPS:-10000}"
MIN_CALLED="${MIN_CALLED:-2}"
THREADS="${THREADS:-4}"
MEMORY_MB="${MEMORY_MB:-8000}"
PYTHON="${PYTHON:-python3}"
PLINK="${PLINK:-plink}"
export AUTOSOME_MAP GROUP_COLUMN OUT_DIR META_PATH FST_DIR BLOCK_BP N_BOOT SEED CHUNK_SNPS MIN_CALLED
# Avoid spawning many BLAS threads for small bootstrap matrix operations.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

command -v "$PYTHON" >/dev/null || { echo "Python 3 is required."; exit 1; }
command -v Rscript >/dev/null || { echo "Rscript with package ape is required."; exit 1; }
"$PYTHON" "${SCRIPT_DIR}/python/fst_check.py"
Rscript --vanilla -e 'if (!requireNamespace("ape", quietly=TRUE)) stop("Install R package ape before running this pipeline.")'
[[ -s "$META_PATH" ]] || { echo "Missing metadata: $META_PATH"; exit 1; }
[[ ! -e "$OUT_DIR" ]] || { echo "OUT_DIR already exists; choose a fresh output directory."; exit 1; }
mkdir -p "$OUT_DIR/plots" "$OUT_DIR/trees"
exec > >(tee "$OUT_DIR/pipeline.log") 2>&1
echo "Output directory: $OUT_DIR"
"$PYTHON" "$SCRIPT_DIR/python/provenance.py" "$SCRIPT_DIR" "$OUT_DIR"
printf 'VCF_PATH=%s\nBFILE=%s\nMETA_PATH=%s\nBLOCK_BP=%s\nN_BOOT=%s\nSEED=%s\nMIN_CALLED=%s\n' \
    "$VCF_PATH" "$BFILE" "$META_PATH" "$BLOCK_BP" "$N_BOOT" "$SEED" "$MIN_CALLED" > "$OUT_DIR/settings.txt"

if [[ -z "$BFILE" ]]; then
    echo "Step 1: Convert the QC VCF to autosomal biallelic SNP genotypes."
    command -v "$PLINK" >/dev/null || { echo "PLINK 1.9 is required for VCF conversion."; exit 1; }
    version=$("$PLINK" --version)
    [[ "$version" == *"PLINK v1.9"* ]] || { echo "PLINK 1.9 required; found: $version"; exit 1; }
    [[ -s "$VCF_PATH" ]] || { echo "Missing VCF: $VCF_PATH"; exit 1; }
    echo "$version"
    chromosomes=$(awk '!/^#/ && NF==2 {print $1; print $2}' "$AUTOSOME_MAP" | paste -sd,)
    "$PLINK" --vcf "$VCF_PATH" --double-id --vcf-half-call missing \
        --chr-set 24 no-xy --allow-extra-chr --chr "$chromosomes" --biallelic-only strict \
        --snps-only just-acgt --keep-allele-order --make-bed \
        --threads "$THREADS" --memory "$MEMORY_MB" --out "$OUT_DIR/qc_autosomal"
    ANALYSIS_BFILE="$OUT_DIR/qc_autosomal"
else
    echo "Step 1: Use supplied PLINK genotypes (must be QC data without LD pruning)."
    ANALYSIS_BFILE="$BFILE"
fi
for suffix in bed bim fam; do
    [[ -s "$ANALYSIS_BFILE.$suffix" ]] || { echo "Missing $ANALYSIS_BFILE.$suffix"; exit 1; }
done
export ANALYSIS_BFILE
echo "ANALYSIS_BFILE=$ANALYSIS_BFILE" >> "$OUT_DIR/settings.txt"

echo "Step 2: Compute pairwise FST components and jointly bootstrap genomic blocks."
"$PYTHON" "${SCRIPT_DIR}/python/fst_analysis.py"

echo "Step 3: Make the NJ tree, bootstrap support, heatmap, and CI plot."
Rscript --vanilla "${SCRIPT_DIR}/r/plot_fst.R"
echo "Completed: $OUT_DIR"
echo "Figures are PNG at 600 dpi. Review METHODS.txt and QC tables before interpretation."
