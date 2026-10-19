need "$BCFTOOLS"; need "$VCFTOOLS"
out="$RUN_DIR/01_qc"; mkdir "$out"
"$BCFTOOLS" view --types snps --min-alleles 2 --max-alleles 2 "$VCF" -Oz -o "$out/01.snp.vcf.gz"
"$BCFTOOLS" filter -i "QUAL>=${MIN_QUAL}" "$out/01.snp.vcf.gz" -Oz -o "$out/02.qual.vcf.gz"
qc_input="$out/02.qual.vcf.gz"
if [[ "$APPLY_DP" == 1 ]]; then
    # Single | evaluates each sample; INFO/DP is not the per-genotype depth.
    "$BCFTOOLS" +setGT "$qc_input" -Oz -o "$out/02.depth.vcf.gz" -- \
        -t q -n . -i "FMT/DP<${MIN_DP} | FMT/DP>${MAX_DP} | FMT/DP=\".\""
    qc_input="$out/02.depth.vcf.gz"
fi
"$VCFTOOLS" --gzvcf "$qc_input" --maf "$MAF" --max-missing "$MISS" \
    --recode --recode-INFO-all --out "$out/03.filtered"
# Recompute tags after genotype masking (avoid stale INFO/AC,AN,AF).
"$BCFTOOLS" +fill-tags "$out/03.filtered.recode.vcf" -Oz -o "$out/03.filtered.recode.vcf.gz" -- -t AC,AN,AF,MAF,NS
"$BCFTOOLS" index -t "$out/03.filtered.recode.vcf.gz"
"$BCFTOOLS" stats "$out/03.filtered.recode.vcf.gz" > "$out/qc.stats"
