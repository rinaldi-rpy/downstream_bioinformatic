need "$BCFTOOLS"; require_plink19; require_file "$QC_VCF"
out="$RUN_DIR/02_ld"; mkdir "$out"
# Only accepted autosomal accessions; numeric mapping avoids PLINK skipping extra chromosomes.
awk '!/^#/ && NF==2 {print $1}' "$AUTOSOME_MAP" > "$out/chromosomes.txt"
chromosomes=$(paste -sd, "$out/chromosomes.txt")
"$BCFTOOLS" view -t "$chromosomes" "$QC_VCF" -Ou | \
    "$BCFTOOLS" annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o "$out/autosomal.vcf.gz"
"$BCFTOOLS" index -t "$out/autosomal.vcf.gz"
"$PLINK" --vcf "$out/autosomal.vcf.gz" --double-id --vcf-half-call missing \
    --allow-extra-chr --biallelic-only strict --snps-only just-acgt --keep-allele-order \
    --threads "$THREADS" --make-bed --out "$out/merged.accessions"
"$PYTHON" "$SCRIPT_DIR/python/remap_bim.py" "$out/merged.accessions.bim" "$AUTOSOME_MAP" "$out/merged.bim"
cp "$out/merged.accessions.bed" "$out/merged.bed"; cp "$out/merged.accessions.fam" "$out/merged.fam"
"$PLINK" --bfile "$out/merged" --chr-set 24 no-xy --indep-pairwise "$LD_WINDOW" "$LD_STEP" "$LD_R2" --threads "$THREADS" --out "$out/merged"
"$PLINK" --bfile "$out/merged" --chr-set 24 no-xy --extract "$out/merged.prune.in" \
    --keep-allele-order --make-bed --out "$out/merged.pruned"
