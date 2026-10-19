need "$ADMIXTURE"; require_file "$LD_PREFIX.bed"
out="$RUN_DIR/04_admixture"; mkdir "$out"
n_samples=$(wc -l < "$LD_PREFIX.fam")
(( K_MAX < n_samples )) || { echo 'K_MAX must be smaller than sample count.' >&2; exit 1; }
printf 'K\tReplicate\tSeed\n' > "$out/seeds.tsv"
for ((k=K_MIN;k<=K_MAX;k++)); do
    for ((rep=1;rep<=N_REPS;rep++)); do
        seed=$((SEED + k*N_REPS + rep))
        run="$out/K${k}_rep${rep}"; mkdir "$run"
        for ext in bed bim fam; do cp "$LD_PREFIX.$ext" "$run/input.$ext"; done
        # Numeric chromosome aliases are a local ADMIXTURE copy only.
        "$PYTHON" "$SCRIPT_DIR/python/numeric_bim.py" "$run/input.bim" "$run/chromosome_map.tsv"
        printf '%s\t%s\t%s\n' "$k" "$rep" "$seed" >> "$out/seeds.tsv"
        (cd "$run"; "$ADMIXTURE" --cv="$CV_FOLDS" -j"$THREADS" -s "$seed" input.bed "$k" > admixture.log 2>&1)
        require_file "$run/input.${k}.Q"; require_file "$run/input.${k}.P"
    done
done
"$PYTHON" "$SCRIPT_DIR/python/admixture_summary.py" "$out" "$K_MIN" "$K_MAX" "$N_REPS"
