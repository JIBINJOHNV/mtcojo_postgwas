# mtcojo_postgwas

An end-to-end **multi-trait conditional GWAS pipeline** built around GCTA mtCOJO,
with optional **PostGWAS Docker harmonisation** and **LDSC genetic correlation** analysis.

The pipeline takes GWAS VCF inputs, validates the manifest and reference resources,
converts each GWAS VCF into the GCTA mtCOJO `.ma` input format, runs GCTA mtCOJO,
optionally runs PostGWAS harmonisation, optionally runs LDSC heritability/genetic
correlation analysis, and writes a timestamped log plus a standalone HTML report
covering inputs, conversion, mtCOJO results, PostGWAS outputs, LDSC outputs,
summary tables, overlap plots, Manhattan/Q-Q plots, and conditional significance shifts.

---

## How It Works

The pipeline runs in sequential stages:

### Stage 1 — Variant ID Format Auto-Detection
Before anything else, the pipeline reads your PLINK BIM file and LD score files to detect the variant ID format:

- Samples up to 2,000 IDs from the BIM file and checks what fraction match `rs\d+`
- Samples the SNP column of the first `*.l2.ldscore.gz` file
- **If BIM and LD score files use different formats → pipeline stops with a clear error**
- The detected format (`rsID` or `CHROM_POS_REF_ALT`) is then automatically used for all
  subsequent VCF conversions — you never need to specify this manually

### Stage 2 — BIM Sanitization (if needed)
Many 1000 Genomes reference panels ship with a non-standard 7-column `.bim` file.
GCTA and PLINK reject this format. The pipeline detects 7-column BIM files, extracts the
correct ID column (based on the format detected in Stage 1), and writes a standard 6-column
BIM alongside symlinked `.bed`/`.fam` files. **The original IDs are never modified.**

### Stage 3 — Single-Pass VCF Extraction
For each GWAS VCF, `bcftools query` extracts all necessary fields in a single streaming pass:

| VCF Field | Meaning | GCTA .ma column |
|-----------|---------|-----------------|
| `%ID` | rsID annotation | `SNP` (if rsID format detected) |
| `%CHROM_%POS_%REF_%ALT` | Constructed ID | `SNP` (if CHRPOS format detected) |
| `[%AF]` | Allele frequency | `freq` |
| `[%ES]` | Effect size (beta) | `b` |
| `[%SE]` | Standard error | `se` |
| `[%LP]` | −log₁₀(p) → p = 10^(−LP) | `p` |
| `[%NEF]` | Effective sample size | `N` |
| `[%SI]` | Imputation INFO score | Passed to PostGWAS |

After extraction, variant IDs are validated against the BIM — if overlap is below 1%,
the pipeline aborts with a diagnostic message showing sample IDs from both sides.

### Stage 4 — GCTA mtCOJO Analysis
Constructs a `.mtcojo.list` mapping each trait to its `.ma` file, then launches GCTA:
1. **Univariate LDSC** — SNP heritability per trait
2. **Bivariate LDSC** — genetic correlation between traits
3. **GSMR** — causal effect of covariate on target using GW-significant index SNPs, after HEIDI-outlier pleiotropic SNP removal
4. **mtCOJO conditioning** — produces per-variant conditional betas (`bC`), standard errors (`bC_se`), and p-values (`bC_pval`) in `.mtcojo.cma`

### Stage 5 — PostGWAS Harmonisation (optional, `--run-postgwas`)
Joins mtCOJO results with in-memory target coordinates, writes the 25-column PostGWAS
manifest (`gwas2vcf_input2.tsv`), and launches the `jibinjv/postgwas:1.4` Docker container
to produce harmonised GRCh37/GRCh38 VCF output.

PostGWAS harmonisation uses the companion workflow:
[JIBINJOHNV/postgwas](https://github.com/JIBINJOHNV/postgwas).

> **Requires:** `--defaults /path/to/harmonisation.yaml` and `--resource-folder /path/to/gwas2vcf/`  
> These are validated **before any analysis starts** — if missing, the pipeline exits immediately.

### Stage 6 — LDSC Heritability and Genetic Correlation (optional, `--run-ldsc`)
Runs LDSC SNP heritability (`--h2`) and pairwise LDSC `--rg` genetic correlation between all input traits **and** the
mtCOJO-conditioned target output. Uses the [CBIIT/ldsc](https://github.com/CBIIT/ldsc)
Python 3 fork, installed locally by `install.sh` — **no Docker required**.

**LDSC is independent of PostGWAS.** You can run LDSC without PostGWAS, PostGWAS without LDSC,
or both together. If running both, pass both flags with their respective required inputs.

LDSC pipeline stages:
1. Convert `.ma` / `.mtcojo.cma` → munge input TSV (EZ = b/se, all fields validated)
2. Run `munge_sumstats.py` in parallel (uses HapMap3 SNP list for variant filtering)
3. Run `ldsc.py --h2` for each trait and conditioned output
4. Run `ldsc.py --rg` for all pairwise trait combinations (parallel, batched)
5. Parse all `.log` files → `ldsc_h2_results.csv` and `ldsc_results.csv`

> **Note on `.mtcojo.cma` conditional columns:** If GCTA's GSMR step had insufficient
> GW-significant SNPs (common for small/regional test datasets), the `bC` columns will be
> `null`. The pipeline automatically falls back to marginal `b/se/p` with a clear warning.

### Stage 7 — Standalone HTML Report
The final report is always generated at the end of a successful pipeline run:

```
<output>/<out>_report.html
```

The report contains the executed command, per-step parameters, generated file inventory,
variant counts, top association tables, Venn/UpSet overlap plots, Manhattan and Q-Q plots,
and LDSC heatmaps when available. LDSC and PostGWAS sections are conditional:

- If `--run-ldsc` is enabled, the HTML report opens LDSC h² and pairwise rg summary sections and embeds the rg heatmap when it can be generated.
- If LDSC CSV/log artifacts already exist under the output directory, regenerated HTML detects them and displays the LDSC summary instead of showing LDSC as skipped.
- If `--run-postgwas` is enabled, the report records the PostGWAS parameters and lists the generated harmonisation artifacts in the output inventory.
- If either optional step is not enabled and no artifacts exist, the report explicitly marks that analysis as skipped.

---

## Installation

### Option 1 — Automatic (Recommended)
Detects your OS (macOS or Linux/Ubuntu), installs Miniforge/mamba if needed,
creates the conda environment with `bcftools` and `gcta64`, then clones and installs
[CBIIT/ldsc](https://github.com/CBIIT/ldsc) scripts locally.

```bash
cd /path/to/mtcojo_postgwas
bash install.sh
conda activate mtcojo_postgwas
```

What `install.sh` does:
1. Detects macOS or Linux, installs Miniforge3 if mamba/conda not found
2. Creates the `mtcojo_postgwas` conda environment (`bcftools`, `gcta`, `python≥3.8`, Python plotting libraries, `r-data.table` for fast R plot-table loading)
3. Installs `mtcojo_postgwas` Python package in editable mode
4. Clones `https://github.com/CBIIT/ldsc` → `tools/ldsc/`
5. Installs `ldsc` Python package via pip
6. Verifies all tools (`bcftools`, `gcta64`, `munge_sumstats.py`, `ldsc.py`, `matplotlib`, `matplotlib-venn`, R plotting helpers)

If you already created the environment before HTML Venn/custom Python plotting and fast R `fread()` support was added, update it with:

```bash
conda install -n mtcojo_postgwas matplotlib matplotlib-venn r-data.table
```

### Option 2 — Manual via mamba/conda
```bash
mamba env create -f environment.yml
conda activate mtcojo_postgwas
conda install -n mtcojo_postgwas matplotlib matplotlib-venn r-data.table
pip install -e .
# Then clone CBIIT/ldsc manually if you need LDSC:
git clone https://github.com/CBIIT/ldsc.git tools/ldsc/
```

---

## Resuming Failed or Interrupted Runs

The pipeline is resumable by default. If expected output files already exist and are non-empty, rerunning the same command reuses them and continues from the first missing or incomplete stage. This is useful after a failure in GCTA, PostGWAS, LDSC, or HTML/report plotting.

Examples:

```bash
# Resume using existing non-empty outputs
mtcojo-postgwas ...same arguments...

# Force all stages to regenerate
mtcojo-postgwas ...same arguments... --force

# Rebuild only the final HTML report when upstream outputs are already present
mtcojo-postgwas ...same arguments... --force-report
```

Resume behavior:
- Existing non-empty `.ma` files are reused.
- Existing mtCOJO `.mtcojo.cma` is reused when the trait list has not changed.
- Existing PostGWAS harmonised outputs are reused.
- Existing LDSC `ldsc_results.csv` and `ldsc_h2_results.csv` are reused.
- Existing HTML is reused only if no upstream stage was rerun; use `--force-report` to refresh it.

---

## Running the Bundled Test

The repository includes a self-contained test dataset (chr1:1Mb–2Mb subset):

```bash
conda activate mtcojo_postgwas
bash run_test.sh
```

**Expected output in `test_output/`:**

| File | Description |
|------|-------------|
| `01_gcta_ma_conversion/SCZ.ma` | Target GCTA summary statistics |
| `01_gcta_ma_conversion/BIP.ma` | Covariate GCTA summary statistics |
| `02_gcta_mtcojo_results/scz_bip_test.mtcojo.cma` | Conditional summary statistics (`bC`, `bC_se`, `bC_pval`) |
| `02_gcta_mtcojo_results/scz_bip_test.badsnps` | SNPs excluded due to allele mismatch |
| `00_manifest_and_logs/scz_bip_test.pipeline.log` | Full timestamped pipeline log |
| `05_plots_and_tables/merged_gwas_summary.tsv` | Wide merged summary statistics used for plots and report tables |
| `05_plots_and_tables/*.png` / `*.jpg` | Manhattan, Q-Q, overlap, and optional LDSC heatmap plots |
| `scz_bip_test_report.html` | Standalone HTML dashboard report |

---

## Bundled Test Data

Located at `mtcojo_postgwas/test_data/`:

| File | Description |
|------|-------------|
| `scz_mini.vcf.gz` | PGC3 SCZ GWAS VCF — chr1:1Mb–2Mb (2,729 SNPs) |
| `bip_mini.vcf.gz` | BIP 2024 EUR GWAS VCF — chr1:1Mb–2Mb (2,260 SNPs) |
| `ref_panel_chr1_mini.*` | 1000G EUR reference panel — 6,243 SNPs, 503 individuals |
| `ld_scores/` | EUR LD scores (chr1 subset) |
| `test_manifest.csv` | Manifest template (paths resolved by `run_test.sh`) |

---

## Input Manifest Format

Create a CSV where **Row 1 = target trait**, **Row 2+ = covariate traits**:

```csv
sample_id,file_path,sample_prevalence,population_prevalence
SCZ,/path/to/PGC3_SCZ_GRCh37.vcf.gz,0.408258,0.01
BIP,/path/to/BIP_GRCh37.vcf.gz,0.25,0.02
```

- `sample_prevalence` / `population_prevalence` — used by GCTA for liability-scale heritability. Set to `NA` for quantitative traits.

---

## Usage Examples

### mtCOJO Only
```bash
mtcojo-postgwas \
  -m manifest.csv \
  -d ./output \
  -o scz_bip_mtcojo \
  -b /path/to/1000G_EUR_reference \
  -l /path/to/eur_w_ld_chr \
  --gwas-thresh 1e-5
```

### mtCOJO + LDSC Heritability and Genetic Correlation
LDSC uses `.ma` files and `.mtcojo.cma` directly — **no PostGWAS needed**.
```bash
mtcojo-postgwas \
  -m manifest.csv \
  -d ./output \
  -o scz_bip_mtcojo \
  -b /path/to/1000G_EUR_reference \
  -l /path/to/eur_w_ld_chr \
  --run-ldsc \
  --ldsc-snp-list /path/to/w_hm3.snplist
```

### mtCOJO + PostGWAS Harmonisation
```bash
mtcojo-postgwas \
  -m manifest.csv \
  -d ./output \
  -o scz_bip_mtcojo \
  -b /path/to/1000G_EUR_reference \
  -l /path/to/eur_w_ld_chr \
  --run-postgwas \
  --defaults /path/to/harmonisation.yaml \
  --resource-folder /path/to/postgwas/gwas2vcf/
```

### mtCOJO + PostGWAS + LDSC (all three)
```bash
mtcojo-postgwas \
  -m manifest.csv \
  -d ./output \
  -o scz_bip_mtcojo \
  -b /path/to/1000G_EUR_reference \
  -l /path/to/eur_w_ld_chr \
  --run-postgwas \
  --defaults /path/to/harmonisation.yaml \
  --resource-folder /path/to/postgwas/gwas2vcf/ \
  --run-ldsc \
  --ldsc-snp-list /path/to/w_hm3.snplist \
  --ldsc-n-parallel 4
```

> **ID format is auto-detected.** The pipeline reads the BIM and LD score files at startup
> and automatically selects rsID or CHROM\_POS\_REF\_ALT mode — no flags required.

> **Genome build and SNP ID consistency are mandatory.** All input files should use the
> same genome build/version. GWAS VCF coordinates and alleles, PLINK `.bed/.bim/.fam`
> reference files, LDSC LD score files, HapMap SNP lists, and PostGWAS resources should
> all be from the same build, for example all GRCh37 or all GRCh38. The SNP ID pattern
> must also match between the PLINK BIM and LDSC files: use rsIDs everywhere or
> `CHROM_POS_REF_ALT`-style IDs everywhere; do not mix ID systems.

---

## Standalone HTML Report

Every successful run writes `<out-dir>/<out>_report.html`. This file is self-contained:
plots are embedded as base64 images, so the report can be opened in a browser or shared
with collaborators without requiring a live Python/R session.

The report is designed for post-GWAS review and reproducibility:

| Section | What it shows |
|---------|---------------|
| Input datasets and trait roles | Target/covariate assignment from the manifest |
| Workflow and command execution | Analysis flow, exact CLI command, and per-step parameter matrix |
| Generated output files | Files detected in `00_manifest_and_logs/`, `01_gcta_ma_conversion/`, `02_gcta_mtcojo_results/`, `03_postgwas_harmonisation/`, `04_ldsc_analysis/`, and `05_plots_and_tables/` |
| Variant and top-hit summaries | Variant counts and top SNP tables from `.ma` and `.mtcojo.cma` files |
| Overlap visualizations | Venn and UpSet plots for all variants and P-value thresholds |
| Multi-trait plots | Stacked Manhattan, overlaid Manhattan, separate Q-Q, and combined Q-Q plots |
| LDSC summary | h² table, pairwise rg table, and rg heatmap when `--run-ldsc` output exists |

Conditional behavior:

- With `--run-ldsc`, the LDSC h² and rg sections are populated from `04_ldsc_analysis/ldsc_h2_results.csv`, `04_ldsc_analysis/ldsc_results.csv`, and LDSC log files.
- Without `--run-ldsc`, the LDSC section is shown as skipped unless previous LDSC artifacts are found in the output directory.
- With `--run-postgwas`, the report includes the PostGWAS configuration in the parameter matrix and lists harmonisation outputs from `03_postgwas_harmonisation/`.
- Without `--run-postgwas`, PostGWAS is marked as skipped in the parameter matrix.

---

## Error Messages

**BIM and LD score ID format mismatch:**
```
ERROR  ✘  ID format mismatch between BIM and LD score files!
           BIM format      : RSID
           LD score format : CHRPOS
           Fix: Use matched BIM + LD score resources.
        ━━━  PIPELINE ABORTED  ━━━
```

**Critically low variant ID overlap:**
```
ERROR  ✘  Critically low variant ID overlap: 0.02%  (11 SNPs)
           VCF format   : RSID
           VCF IDs (5)  : ['1_752566_A_G', '1_752721_T_C', ...]
           BIM IDs (5)  : ['rs3094315', 'rs3131972', ...]
        ━━━  PIPELINE ABORTED  ━━━
```

**PostGWAS inputs missing (detected before analysis starts):**
```
ERROR: --run-postgwas is enabled but required inputs are missing:

  --defaults / harmonisation.yaml not found: /path/...
    Pass: --defaults /path/to/harmonisation.yaml

  --resource-folder (gwas2vcf resources) not found: /path/...
    Pass: --resource-folder /path/to/postgwas/gwas2vcf/
```

**LDSC scripts not found (run install.sh first):**
```
ERROR  ✘  LDSC scripts not found.
           Run:  bash install.sh
           Or pass:  --ldsc-dir /path/to/CBIIT_ldsc_clone
```

---

## Command Line Reference

### Core
| Flag | Description | Default |
|------|-------------|---------|
| `-m`, `--manifest` | Input manifest CSV (*required*) | — |
| `-o`, `--out` | Output file prefix | `scz_bip_mtcojo` |
| `-d`, `--out-dir` | Output directory | `./` |

### GCTA mtCOJO
| Flag | Description | Default |
|------|-------------|---------|
| `-b`, `--bfile` | PLINK reference panel prefix | — |
| `--mbfile` | File listing multiple PLINK panels | — |
| `-l`, `--ref-ld-chr` | LD score files directory | — |
| `-w`, `--w-ld-chr` | LD score weights directory | same as `-l` |
| `--gwas-thresh` | Index SNP p-value threshold | `1e-5` |
| `--clump-r2` | LD r² threshold for index SNP clumping | — |
| `--heidi-thresh` | HEIDI pleiotropic SNP removal cutoff | `0.01` |
| `--gsmr-snp-min` | Min index SNPs for GSMR | — |
| `--gsmr-ld-fdr` | FDR threshold for LD pruning in GSMR | — |
| `--diff-freq` | Max allele frequency difference vs BIM | — |
| `--mtcojo-bxy` | Pre-computed causal effect file | — |
| `--no-bim-fix` | Skip 7-column BIM sanitization | `False` |

### Tool Paths
| Flag | Description |
|------|-------------|
| `--bcftools` | Path to `bcftools` binary |
| `--gcta64` | Path to `gcta64` binary |

### PostGWAS Harmonisation (`--run-postgwas`)
> All inputs are validated before any analysis starts.

| Flag | Description | Default |
|------|-------------|---------|
| `--run-postgwas` | Enable PostGWAS Docker harmonisation step | `False` |
| `--defaults` | Path to `harmonisation.yaml` (*required*) | — |
| `--resource-folder` | gwas2vcf reference resources directory (*required*) | — |
| `--sdd-dir` | Host directory mounted into Docker | auto |
| `--liftover` | GRCh37→GRCh38 liftover (`Yes`/`No`) | `No` |
| `--docker-image` | Docker image | `jibinjv/postgwas:1.4` |
| `--docker-platform` | Docker platform passed to `docker run` | `linux/amd64` |
| `--nthreads` | Docker thread count | `23` |
| `--max-mem` | Docker max memory | `50G` |

### LDSC Heritability and Genetic Correlation (`--run-ldsc`)
> Independent of `--run-postgwas`. Uses `.ma` files and `.mtcojo.cma` directly.

| Flag | Description | Default |
|------|-------------|---------|
| `--run-ldsc` | Enable LDSC pairwise rg analysis | `False` |
| `--ldsc-snp-list` | HapMap3 SNP list (w_hm3.snplist) (*required*) | — |
| `--ldsc-ld-dir` | LD reference directory | same as `-l` |
| `--ldsc-dir` | Path to CBIIT/ldsc clone (auto-detected after `install.sh`) | auto |
| `--ldsc-n-parallel` | Parallel processes for munge/rg | `4` |
| `--ldsc-batch-size` | Target traits per `ldsc.py --rg` call | `10` |
| `--ldsc-center-z` | Center Z-scores to median 0 for small/regional datasets that fail munge checks | `False` |

---

## Output Directory Layout

The pipeline uses stable numbered subdirectories so logs, intermediate files, optional
analyses, and report assets remain easy to inspect:

```
output/
├── 00_manifest_and_logs/
│   └── <out>.pipeline.log
├── 01_gcta_ma_conversion/
│   ├── <trait>.ma
│   └── <trait>_commands.sh
├── 02_gcta_mtcojo_results/
│   ├── <out>.mtcojo.list
│   ├── <out>.mtcojo.cma
│   ├── <out>.badsnps
│   ├── <out>.log
│   └── <out>_commands.sh
├── 03_postgwas_harmonisation/       # only populated when --run-postgwas succeeds
├── 04_ldsc_analysis/                # only populated when --run-ldsc succeeds
├── 05_plots_and_tables/
│   ├── merged_gwas_summary.tsv
│   ├── multi_trait_manhattan.png
│   ├── multi_trait_qq.png
│   ├── multi_trait_qq_combined.png
│   └── overlap / optional R plot files
└── <out>_report.html
```

## LDSC Output Files

When `--run-ldsc` is enabled, the following are written under the output directory:

```
output/
└── 04_ldsc_analysis/
    ├── ldsc_munge_input/       {trait}_mungeinput.tsv  (one per trait + conditioned output)
    ├── ldsc_sumstats/          {trait}.sumstats.gz      (after munge_sumstats.py)
    ├── ldsc_h2/                {trait}.log              (ldsc.py --h2 output)
    ├── ldsc_results/           {ref}_batch*.log         (ldsc.py --rg output)
    ├── ldsc_h2_results.csv     per-trait SNP heritability summary
    ├── ldsc_results.csv        all pairwise rg values compiled
    └── ldsc_rg_heatmap.png     rg heatmap when plotting succeeds
```

**`ldsc_h2_results.csv` columns:**

| Column | Description |
|--------|-------------|
| `trait` | Trait name |
| `n_snps` | Number of SNPs used by LDSC |
| `h2_observed` | Observed-scale SNP heritability |
| `h2_observed_se` | Standard error for observed-scale h² |
| `h2_liability` | Liability-scale h² when prevalence values are available |
| `h2_liability_se` | Standard error for liability-scale h² |
| `intercept` | LDSC regression intercept |
| `lambda_gc` | Genomic control inflation factor |
| `mean_chi2` | Mean chi-square statistic |

**`ldsc_results.csv` columns:**

| Column | Description |
|--------|-------------|
| `p1`, `p2` | Trait names |
| `rg` | Genetic correlation estimate |
| `se` | Standard error |
| `z` | Z-score |
| `p` | P-value |
| `h2_obs` | Observed-scale heritability |
| `gcov_int` | Genetic covariance intercept |

## PostGWAS Output Files

When `--run-postgwas` is enabled, outputs are written under:

```
output/
└── 03_postgwas_harmonisation/
    ├── gwas2vcf_input2.tsv        25-column PostGWAS harmonisation input
    ├── config.yaml                generated PostGWAS run manifest
    └── sdd_disk/                  Docker-mounted working directory unless --sdd-dir is provided
```

The exact files produced by the Docker container depend on the PostGWAS resource bundle
and harmonisation defaults. The HTML report inventories whatever files are present after
the run completes.

---

## Development and GitHub Workflow

This repository is easiest to work on from the conda environment used by the pipeline:

```bash
conda activate mtcojo_postgwas
pip install -e .
python -m compileall -q mtcojo_postgwas
bash run_test.sh
```

Recommended contribution workflow:

1. Open an issue or describe the scientific/software change clearly before broad refactors.
2. Keep pull requests focused: one bug fix, feature, or documentation update per PR.
3. Include the command used for validation in the PR description.
4. For scientific logic changes, document the affected stage, expected output files, and any assumptions about GCTA, LDSC, PostGWAS, build, or reference resources.
5. Do not commit local output directories, large reference resources, `.DS_Store`, conda environments, or generated logs.

Before opening a PR, run at least:

```bash
python -m compileall -q mtcojo_postgwas
bash run_test.sh
```

---

## Package Structure

```
mtcojo_postgwas/
├── install.sh              ← Auto-installer (macOS + Linux, clones CBIIT/ldsc)
├── run_test.sh             ← One-command bundled test
├── environment.yml         ← Conda env: bcftools, gcta, python deps, ldsc
├── pyproject.toml
├── README.md
├── tools/                  ← Created by install.sh
│   ├── ldsc/               ← CBIIT/ldsc clone (munge_sumstats.py, ldsc.py)
│   └── ldsc_paths.txt      ← Script paths written by install.sh
└── mtcojo_postgwas/
    ├── __init__.py
    ├── cli.py              ← CLI entry point, arg parsing, stage orchestration
    ├── core/
    │   └── logger.py       ← Colourful stdout + file logging, ASCII logo
    ├── io/
    │   ├── bim_sanitizer.py
    │   │                   ← BIM ID detection, LD score cross-check, 7→6-col fix
    │   └── vcf_converter.py
    │                       ← Single-pass bcftools extraction, overlap validation
    ├── stages/
    │   ├── gcta.py         ← GCTA mtCOJO subprocess wrapper
    │   ├── postgwas.py     ← PostGWAS manifest builder + Docker launcher
    │   └── ldsc.py         ← LDSC munge + h²/rg pipeline (CBIIT/ldsc, no Docker)
    ├── reporting/
    │   ├── report_generator.py
    │   │                   ← Standalone HTML report, tables, and plot orchestration
    │   ├── conditional_shift_summary.py
    │   │                   ← Gained/lost significance summaries after conditioning
    │   └── assets/
    │       └── plot_gwas.R ← Optional rMVP/CMplot plotting helper
    └── test_data/
        ├── scz_mini.vcf.gz
        ├── bip_mini.vcf.gz
        ├── ref_panel_chr1_mini.*
        ├── ld_scores/
        └── test_manifest.csv
```

---

## License

MIT License
