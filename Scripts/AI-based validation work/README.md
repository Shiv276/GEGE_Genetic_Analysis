# GEGE3888 — Exploratory Bd phenotype GWAS

This is **supplementary / exploratory analysis**, not an independent replication or proof of resistance loci. The scripts are separate from the established infection-establishment GWAS and **do not overwrite it**.

## What it tests

| Analysis | Phenotype | Approx. n | Statistical model |
|---|---|---:|---|
| `survival98` | `survived_98` = 1 if no recorded euthanasia endpoint through day 98, else 0 | 892 (359 survived, 533 endpoint) | Binomial logistic GMMAT with GRM |
| `early_detection` | `early_detected` = 1 if **any** Bd ITS measurement >0 in weeks 1–6 | 880 (566 positive, 314 negative) | Binomial logistic GMMAT with GRM |
| `early_positive_burden` | `log10(maximum weeks 1–6 Bd zoospore equivalents + 1)` **only among early-positive frogs** | 566 | Gaussian GMMAT with GRM; diagnostics necessary |

The first analysis is a **binary endpoint** comparison, not a time-to-event GWAS. Its use of the observed `euthanised` indicator assumes it represents the clinically relevant endpoint as recorded in the published experiment. An eventual survival-time analysis should account explicitly for censoring and family relatedness.

The early outcome uses the six-week window because 880/892 frogs have all six weekly results; no missing value is treated as a zero. The conversion factor 86 ITS target copies per zoospore equivalent and Wastepoint strain come from `Psco_RepeatedExp.R`. **Any positive** is a qPCR detection phenotype; it is *not* the authors' at-least-two-positives infection-establishment definition. The conditionally positive burden phenotype is vulnerable to selection/collider bias, and log-transforming does not guarantee Gaussian residuals. Use these as hypothesis-generating sensitivity analyses, not proof of causality.

## Folder setup

Place the **three scripts and this README** in exactly:

```text
Genetic Data/
└── Scripts/
    ├── config.py
    └── AI-based validation work/
        ├── 01_prepare_phenotypes.py
        ├── 02_run_gwas.R
        ├── 03_compare_results.py
        └── README.md
```

The Python scripts use the existing `Scripts/config.py` value `PHENOTYPE_DIR`. They determine `GWAS_ROOT = PHENOTYPE_DIR.parent / 'gwas'`. The following **already generated** project files must exist:

```text
data/processed/phenotypes/matched_phenotypes.csv       # Includes raw Bd_wk1...Bd_wk6

data/processed/gwas/gwas_input_with_pcs.tsv            # From your earlier four-model R script

data/processed/gwas/grm_892.csv

data/processed/gwas/grm_892_ids.tsv

data/processed/gwas/infection_qc.bed

data/processed/gwas/infection_qc.bim

data/processed/gwas/infection_qc.fam
```

**Do not use** `bd_phenotypes_preliminary.csv` in place of `matched_phenotypes.csv`: the preliminary table discards the individual weekly Bd measurements.

All new results are created under:

```text
data/processed/gwas/ai_validation_work/
    analysis_manifest.tsv
    null_model_summary.tsv
    all_gwas_results_summary.csv
    survival98/
    early_detection/
    early_positive_burden/
```

Each outcome folder contains its own PLINK QC dataset, phenotype table, sample-aligned subset GRM and association/figure outputs. Because the **sample sizes differ between outcomes**, SNP MAC/missingness QC is recalculated for each outcome; SNP sets may therefore differ *between* outcomes, but must be identical *between models of the same outcome*.

## Installation

Python environment you already use:

```powershell
python -m pip install numpy pandas scipy matplotlib statsmodels
```

PLINK2 must be runnable by `plink2` or defined in `Scripts/config.py` as `PLINK2` (the executable path). In RStudio run **once** if needed:

```r
install.packages("GMMAT")
```

No additional Python GWAS software is necessary: GMMAT does the mixed-model association testing in R.

## Exact run order

**1. In VS Code / PowerShell**, open the folder containing the new scripts:

```powershell
cd "C:\Users\Shive\OneDrive\Desktop\GEGE3888\Genetic Data\Scripts\AI-based validation work"
python ".\01_prepare_phenotypes.py"
```

This matches and checks IDs, constructs the outcomes, subsets GRMs, and runs PLINK2 separately for the three sample sets with `--geno 0.05 --mac 20`. The original 39,369-marker infection QC data are the starting point, not the LD-pruned kinship panel. The script prints final retained SNP counts (which cannot be known until PLINK has run on each subset).

For phenotype/GRM checks *without* running PLINK2, use:

```powershell
python ".\01_prepare_phenotypes.py" --skip-plink
```

The R step still requires the PLINK results, so omit `--skip-plink` for the real analysis.

**2. In RStudio**, open `02_run_gwas.R`, click **Source**, then select the generated **`data/processed/gwas/ai_validation_work/analysis_manifest.tsv`** in the file dialog. You do not need to change `getwd()` or `setwd()`.

Default models for each trait:

- GRM + experiment (recommended minimally adjusted model)
- GRM + experiment + PC1–PC3 (PCA sensitivity model)

Set `RUN_ALL_FOUR <- TRUE` in the R script if you also want GRM-only and GRM+PC-only models, as in your original GWAS. All results include the same correctly ID-aligned GRM random effect. Convergence and resulting variant count are checked; fitted models are saved as RDS files for audit.

**3. Back in VS Code / PowerShell**:

```powershell
python ".\03_compare_results.py"
```

Outputs include per-outcome QQ overlays, Manhattan figures, phenotype histograms, top-SNP comparisons, per-model BH-FDR and Bonferroni adjustments, and the aggregate `all_gwas_results_summary.csv`.

## Scientific guardrails for reporting

1. **Do not pick significant-looking models after seeing the results and call them confirmatory.** These are three exploratory phenotypes of overlapping animals, with multiple SNPs and multiple models. Report the full set of tested outcomes and model specifications. Within-model BH corrections are not a correction across all outcomes and model choices.
2. **No association passing FDR or Bonferroni is a valid result.** Smallest nominal p-value and Manhattan peaks are not biological validation, particularly in strongly related families.
3. **A λGC closer to 1 is not sufficient evidence of a superior model.** Inspect QQ-tail behaviour, convergence, covariate justification, and signal stability. Experiment is included because experimental infection rates differ markedly; three PCs are an exploratory sensitivity adjustment.
4. **Early detections are not clinical resistance:** Any-positive by week 6 is not the original two-positive criterion and six weeks are not necessarily independent of events in that period.
5. **Conditioning on early-positive frogs selects a nonrandom subgroup.** Burden associations in the 566 frogs are conditional and may be biased. Check distribution/residual assumptions and do not treat them as unconditional burden effects.
6. The genomic matrix is the **previously computed LD-pruned GRM**, subsetted to each phenotype's sample. GWAS tests use the broader filtered SNP set. Same-chromosome SNPs in the GRM can produce proximal contamination (LOCO GRMs not implemented in this screen).
7. If an association seems promising, it needs independent validation, attention to effect size, possible family/experiment confounding, and ideally a better-powered external cohort before it could inform marker-assisted conservation breeding.

## Interpreting the deliverables

Each phenotype folder generates:

- `phenotype.tsv`, `keep.txt`, `grm.csv`, `grm_ids.tsv`: exactly aligned sample records and relatedness matrix.
- `snps_qc.bed/.bim/.fam/.log`: QC-subset genotype panel and PLINK log.
- `null_<model>.rds`: fitted null GLMM/LMM.
- `gwas_<model>.tsv`: GMMAT per-SNP score-test association table.
- `processed_<model>.csv`: p-values plus within-scan FDR/Bonferroni adjustment.
- `qq_models.png/.pdf`, `manhattan_models.png/.pdf`, `phenotype_distribution.png`, `top_snps_across_models.csv`.

Inspect **`all_gwas_results_summary.csv` first**. Keep phenotype definitions, sample sizes and QC filtering in the report even if no SNP survives correction.

## Data provenance

The original phenotype workbook `Overall-infection.xlsx` and `Psco_RepeatedExp.R` are the source of the challenge data and the Wastepoint ITS-to-zoospore conversion. The existing matched frog data and GRM are the outputs of your earlier, explicitly documented GEGE3888 pipeline. The early-window and conditional-positive outcomes are **new exploratory derived phenotypes**, not presented as the published authors' phenotypes.
