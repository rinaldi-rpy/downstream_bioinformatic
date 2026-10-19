# Buffalo population genomics workflow

Code revision: 19 October 2026.

This workflow covers SNP QC, LD pruning, PCA, replicated ADMIXTURE,
windowed pairwise Weir–Cockerham FST, shared genomic-block FST bootstrap
with an optional NJ visualization, and replicated TreeMix models.
The original scripts, including commented stages and alternative analyses,
are preserved byte-for-byte under `legacy/original/` with a SHA-256 manifest.
A commented stage in the originals is **not** evidence that it was unused.

## Run all analyses with one command

```bash
bash run_all.sh /path/Asia_buff_all.vcf.gz /path/metadata.csv
```

`run_all.sh` runs QC → LD pruning → PCA → replicated ADMIXTURE → pairwise
FST → FST block bootstrap/NJ → replicated TreeMix and its summaries/plots.
Keep the entire workflow folder together: Python/R remain separate files,
but no manual invocation of those files is required. Tool executables, NumPy,
and R packages ape/RColorBrewer are checked before QC starts. The launcher
uses the current revised parameters without changing filtering or algorithms.
Each successful stage feeds the next; an error stops the pipeline for inspection.
Use `bash run_all.sh --help` for usage. It does not install dependencies,
submit a cluster scheduler job, or automatically choose TreeMix migration count.
The production-validation limitations below still apply.

Final folder: `/home/rizalrinaldi/Documents/downstream_github`. See [file naming](docs/file_naming.md) for the old-to-new name mapping.

## Layout

- Root `.sh`: public entry points and stage selection.
- `shell/`: Bash stages and shared configuration loading.
- `python/`: metadata validation, genotype analysis, conversion and summaries.
- `r/`: plotting only, including the supplied local TreeMix helper with documented repairs.
- `config/parameters.sh`: resolved defaults; environment variables can override them.
- `config/autosomes.tsv`: explicit accession-to-numeric chromosome mapping.
- `tests/`: synthetic regression and orchestration tests.
- `docs/audit.md`: file-by-file audit, behavior changes and scientific limitations.
- `docs/validation.md`: what was actually tested and what still needs server validation.
- `legacy/`: historical reference material; never called by the revised workflow.

No R or Python analysis is embedded in the active shell scripts. No active script
installs packages or downloads executable code during analysis.

## Dependencies

Use Bash 4+, Python 3.9+ with NumPy >=1.20, bcftools with `setGT` and `fill-tags`
plugins, VCFtools, **PLINK 1.9**, ADMIXTURE and TreeMix. R requires `ape` for
FST/NJ plots and `RColorBrewer` for TreeMix plots. Excel metadata additionally
requires pandas and an appropriate Excel reader; CSV/TSV uses the standard library.
Install dependencies separately using your cluster's supported environment.
Record exact versions in the environment used for the final results.
The revised LD stage uses PLINK 1.9 consistently; the original used PLINK 2 for LD
and PLINK 1.9 elsewhere. Results should be compared before replacing published outputs.

## Run from the code directory on the cluster

```bash
bash run_pipeline.sh --help

# Default stages: qc,ld,pca,admixture,fst
bash run_pipeline.sh /path/Asia_buff_all.vcf.gz metadata.csv

# All requested analyses (including FST bootstrap and TreeMix)
STEPS=qc,ld,pca,admixture,fst,bootstrap,treemix \
  bash run_pipeline.sh /path/Asia_buff_all.vcf.gz metadata.csv
```

A fresh directory under `runs/analysis_<timestamp>_<pid>/` holds the outputs.
Set `RUN_DIR=/absolute/new/directory` to choose it. Existing output directories
are refused. Stage outputs include `01_qc`, `02_ld`, `03_pca`, `04_admixture`,
`05_fst`, `05_bootstrap`, and `06_treemix` when selected.
Run stages in dependency order; omitted prerequisites must exist under `INPUT_RUN`.

```bash
# Repeat selected analyses without commenting code or overwriting earlier outputs.
INPUT_RUN=/path/to/previous/analysis_run STEPS=pca,admixture \
  bash run_pipeline.sh /path/Asia_buff_all.vcf.gz metadata.csv

# Repeat FST + bootstrap using the unpruned 02_ld/autosomal.vcf.gz.
INPUT_RUN=/path/to/previous/analysis_run STEPS=fst,bootstrap \
  bash run_pipeline.sh /path/Asia_buff_all.vcf.gz metadata.csv

# Reuse OLD unpruned QC VCF directly for bootstrap, independent of the new layout.
VCF_PATH=/path/work/01_qc/03.filtered.recode.vcf.gz \
META_PATH=/path/metadata.csv N_BOOT=2000 BLOCK_BP=1000000 \
  bash run_fst_bootstrap.sh

# Use the PLINK prefix already used in the original TreeMix analysis.
BFILE=/path/work/02_ld/merged.pruned META_PATH=/path/metadata.csv \
MAX_M=8 REPS=10 BLOCK_SIZE=500 bash run_treemix.sh

# Only repeat summary/plots of an existing replicated TreeMix run.
bash summarize_treemix.sh /path/to/treemix/run_directory

# Explicit QC summary; refuses an existing output CSV.
bash summarize_qc.sh /path/raw.vcf.gz /path/qc.vcf.gz \
  /path/merged.pruned.bim /path/new_summary.csv
```

Legacy `INPUT_RUN` directories must expose the documented stage filenames.
For old layouts that lack `02_ld/autosomal.vcf.gz`, run the revised `ld` stage
from the old QC directory or use the independent bootstrap entry point.
The resume entry point expects `runs/treemix_m<M>_rep<R>.llik` and `.treeout.gz`;
it does not infer replicate stability from the old single-model outputs.

## Parameters and deliberate behavior changes

| Setting | Default | Meaning |
|---|---:|---|
| THREADS | 32 | Tool thread limit; independent runs are sequential |
| MIN_QUAL | 20 | Minimum site QUAL |
| MIN_DP / MAX_DP | 5 / 100 | Per-genotype FORMAT/DP limits, inclusive |
| APPLY_DP | 1 | Mask low/high/missing-DP genotypes before site filtering |
| MAF | 0.05 | Global minor allele frequency cutoff |
| MISS | 0.90 | Minimum site call rate, i.e. maximum 10% missing |
| LD_WINDOW / LD_STEP / LD_R2 | 50 / 10 / 0.2 | Window in SNPs, step in SNPs, r-squared cutoff |
| N_PC | 253 | Requested count; effective count <= min(samples-1, SNPs), recorded |
| K_MIN / K_MAX / N_REPS | 2 / 10 / 10 | ADMIXTURE K range and independent starts per K |
| CV_FOLDS | 5 | ADMIXTURE cross-validation folds |
| FST_WINDOW / FST_STEP | 50000 / 10000 | Window/step in bp for descriptive windowed FST |
| N_BOOT / BLOCK_BP | 2000 / 1000000 | FST bootstrap replicates / physical block length |
| MIN_CALLED | 2 | Complete diploid individuals per population at each FST SNP |
| MAX_M / REPS / BLOCK_SIZE | 8 / 10 / 500 | TreeMix migration range 0..8 / starts per m / SNP block size |
| SEED | 20261015 | Reproducible base seed; derived run seeds are saved |
| ROOT | empty | TreeMix outgroup; set only if scientifically justified |
| GROUP_COLUMN | Group Label | Metadata column defining populations |

The dedicated original TreeMix pipeline used `MAX_M=8`; the original main script
used 5. Choose `MAX_M=5` explicitly to retain the latter range. TreeMix uses the
LD-pruned prefix by default to preserve the dedicated original input choice;
set `BFILE` explicitly for an unpruned sensitivity analysis. Keep this choice
and `BLOCK_SIZE` in the reported methods.

The old QC never used its DP settings. `APPLY_DP=0` preserves that omission for
a comparison run; `APPLY_DP=1` is a scientific filtering change and requires QC
and dependent analyses to be rerun. Missing FORMAT/DP is not substituted with
INFO/DP. INFO allele-count tags are recomputed after filtering.
No new HWE, individual call-rate, relatedness or genotype-quality threshold is
introduced without an existing parameter. QC here starts from a merged VCF;
it does not implement raw-read QC, alignment or variant calling.

The default chromosome map retains `NC_059157.1`–`NC_059180.1`, as in the supplied
FST bootstrap and dedicated TreeMix scripts. **It was inherited from the code,
not confirmed against a supplied reference assembly**. Inspect your reference
assembly/VCF header before production use; change `AUTOSOME_MAP` if necessary.
QC retains all biallelic SNP chromosomes; downstream LD/PCA/ADMIXTURE/FST/TreeMix
use the selected autosomes. Numeric PLINK aliases 1–24 are explicit and never
assigned by arbitrary encounter order for LD. X/MT are excluded downstream.
The mapping currently supports up to 24 autosomes.

## Statistical interpretation and reproducibility

FST bootstrap recomputes diploid Weir–Cockerham components from **unpruned**
genotypes and uses `sum(a)/sum(a+b+c)`. Fixed physical blocks do not cross
chromosomes. Blocks with selected SNPs are resampled within chromosome, with
the same draws for all population pairs. Confidence intervals are percentile
95%; negative FST values remain in tables and CIs. Overlapping FST windows
are not independent bootstrap units. Missing-call requirements can cause
differences from VCFtools window/global log estimates; compare site sets.
Bootstrap conditions on the sampled individuals and ascertainment filters.
One Mb is an inherited starting value; report block-length sensitivity.

The NJ tree is an unrooted descriptive summary of FST. Negative input distances
are truncated only when constructing NJ trees; tables retain them. Support
measures stability of splits under genomic-block resampling, not migration
support or a species phylogeny. Four populations are required by the NJ workflow.

ADMIXTURE writes each K/replicate to a separate directory with a saved seed,
then keeps CV errors with replicate identities and selects minimum CV per K.
Inspect distributions and convergence before interpreting K. Q columns can
permute across starts; Q matrices are **not averaged** without label alignment.
No inference of ancestry proportions from unaligned replicate averages is made.

TreeMix preserves genomic SNP order, uses exact A1/A2 counts with consistent
allele orientation, drops SNPs without observed alleles in every population,
and retains population label mappings. Repeated seeds assess optimization
stability; they are **not bootstrap replicates**. No automatic migration count
or migration-direction claim is produced. Likelihood, replicate variation,
residuals and biological rooting must be reviewed. A TreeMix block of 500 SNPs
is different from a one-Mb FST block; the TreeMix format has no chromosome
boundary field, so fixed SNP blocks may straddle concatenated boundaries.

Each new analysis/model run records parameters, code hashes, available tool
versions and paths in `provenance.json`; per-tool logs and seed tables are kept.
Freeze the final environment and add original VCF/assembly accession/checksum,
software versions for ADMIXTURE/TreeMix, and final output checksums to the
submission record. Source hashes alone do not establish input-data identity.

## Validation

```bash
python3 -m unittest discover -s tests -v
Rscript --vanilla tests/test_plots.R
```

See `docs/validation.md` for exact limits. These tests do not replace a full
cluster rerun using the real VCF and installed bioinformatics binaries.
