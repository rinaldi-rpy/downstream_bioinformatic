need "$VCFTOOLS"; need "$BCFTOOLS"; require_file "$AUTO_VCF"
out="$RUN_DIR/05_fst"; mkdir "$out"
"$BCFTOOLS" query -l "$AUTO_VCF" > "$out/samples.txt"
"$PYTHON" "$SCRIPT_DIR/python/metadata.py" populations --meta "$META" --group "$GROUP_COLUMN" --samples "$out/samples.txt" --out "$out/pop"
shopt -s nullglob
popfiles=("$out"/pop/*.txt)
(( ${#popfiles[@]} >= 2 )) || { echo 'FST requires at least two populations.' >&2; exit 1; }
for ((i=0;i<${#popfiles[@]}-1;i++)); do
    for ((j=i+1;j<${#popfiles[@]};j++)); do
        pop1=$(basename "${popfiles[i]}" .txt); pop2=$(basename "${popfiles[j]}" .txt)
        "$VCFTOOLS" --gzvcf "$AUTO_VCF" --weir-fst-pop "${popfiles[i]}" --weir-fst-pop "${popfiles[j]}" \
            --fst-window-size "$FST_WINDOW" --fst-window-step "$FST_STEP" --out "$out/${pop1}_vs_${pop2}"
    done
done
