require_plink19; require_file "$LD_PREFIX.bed"
out="$RUN_DIR/03_pca"; mkdir "$out"
n_samples=$(wc -l < "$LD_PREFIX.fam"); n_variants=$(wc -l < "$LD_PREFIX.bim")
effective_pc=$N_PC
(( effective_pc < n_samples )) || effective_pc=$((n_samples-1))
(( effective_pc <= n_variants )) || effective_pc=$n_variants
(( effective_pc >= 1 )) || { echo 'PCA needs at least two samples and one SNP.' >&2; exit 1; }
printf 'Requested_PCs\tEffective_PCs\tSamples\tSNPs\n%s\t%s\t%s\t%s\n' "$N_PC" "$effective_pc" "$n_samples" "$n_variants" > "$out/pca_dimensions.tsv"
if (( effective_pc != N_PC )); then echo "PCA: requested $N_PC PCs, using $effective_pc for $n_samples samples/$n_variants SNPs."; fi
"$PLINK" --bfile "$LD_PREFIX" --chr-set 24 no-xy --allow-extra-chr --pca "$effective_pc" --threads "$THREADS" --out "$out/pca"
