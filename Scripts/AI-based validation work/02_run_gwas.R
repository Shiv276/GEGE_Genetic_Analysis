# GEGE3888 | Exploratory Bd survival and early-burden GWAS
# Run after: python 01_prepare_phenotypes.py
# Requires the R package GMMAT (install.packages("GMMAT") once).
# Select 'analysis_manifest.tsv' in the dialog. R working directory is irrelevant.

library(GMMAT)

root <- dirname(file.choose())  # Select analysis_manifest.tsv in ai_validation_work/
manifest <- read.delim(file.path(root, "analysis_manifest.tsv"),
                       stringsAsFactors = FALSE)
required <- c("analysis", "trait", "family", "n_frogs")
stopifnot(all(required %in% names(manifest)))

# Default: two experiment-adjusted models per phenotype.
# Optional: set RUN_ALL_FOUR <- TRUE to also fit intercept-only and PC-only.
RUN_ALL_FOUR <- FALSE

summaries <- list()
for (i in seq_len(nrow(manifest))) {
  name <- manifest$analysis[i]
  trait <- manifest$trait[i]
  likelihood <- manifest$family[i]
  folder <- file.path(root, name)
  cat("\n=======================================\nANALYSIS:", name, "\n")

  pheno <- read.delim(file.path(folder, "phenotype.tsv"),
                      stringsAsFactors = FALSE)
  ids <- read.delim(file.path(folder, "grm_ids.tsv"),
                    stringsAsFactors = FALSE)
  K <- as.matrix(read.csv(file.path(folder, "grm.csv"),
                          header = FALSE, check.names = FALSE))
  storage.mode(K) <- "double"

  # Important: compare numeric GRM values; R-generated V1/V2 names differ.
  stopifnot(
    nrow(pheno) == manifest$n_frogs[i],
    !anyDuplicated(pheno$IID),
    identical(as.character(pheno$IID), as.character(ids$IID)),
    all(dim(K) == c(nrow(pheno), nrow(pheno))),
    all(is.finite(K)),
    max(abs(K - t(K))) < 1e-8,
    !anyNA(pheno[[trait]]),
    all(is.finite(as.matrix(pheno[c("PC1", "PC2", "PC3")]))),
    file.exists(file.path(folder, "snps_qc.bed")),
    file.exists(file.path(folder, "snps_qc.bim")),
    file.exists(file.path(folder, "snps_qc.fam"))
  )
  if (likelihood == "binomial") {
    stopifnot(all(pheno[[trait]] %in% c(0, 1)),
              length(unique(pheno[[trait]])) == 2L)
    family_arg <- binomial(link = "logit")
  } else if (likelihood == "gaussian") {
    stopifnot(is.numeric(pheno[[trait]]), all(is.finite(pheno[[trait]])),
              stats::sd(pheno[[trait]]) > 0)
    family_arg <- gaussian(link = "identity")
  } else {
    stop("Unexpected outcome family: ", likelihood)
  }

  fam <- read.table(file.path(folder, "snps_qc.fam"), header = FALSE,
                    stringsAsFactors = FALSE)
  stopifnot(nrow(fam) == nrow(pheno),
            setequal(as.character(fam$V2), as.character(pheno$IID)))

  rownames(K) <- pheno$IID
  colnames(K) <- pheno$IID
  pheno$experiment <- factor(pheno$experiment)
  stopifnot(nlevels(pheno$experiment) == 3L)

  # A and C provided for a direct GRM vs GRM+PC comparison if desired.
  models <- list(
    grm_only = as.formula(paste(trait, "~ 1")),
    grm_experiment = as.formula(paste(trait, "~ experiment")),
    grm_pc3 = as.formula(paste(trait, "~ PC1 + PC2 + PC3")),
    grm_experiment_pc3 = as.formula(
      paste(trait, "~ experiment + PC1 + PC2 + PC3"))
  )
  to_run <- if (RUN_ALL_FOUR) names(models) else
    c("grm_experiment", "grm_experiment_pc3")

  for (model_name in to_run) {
    cat("\n", name, "/", model_name, "\n")
    fit <- glmmkin(
      fixed = models[[model_name]], data = pheno,
      kins = K, id = "IID", family = family_arg
    )
    if (!isTRUE(fit$converged)) {
      stop("Null model failed to converge: ", name, " / ", model_name)
    }
    saveRDS(fit, file.path(folder, paste0("null_", model_name, ".rds")))

    outfile <- file.path(folder, paste0("gwas_", model_name, ".tsv"))
    glmm.score(
      obj = fit,
      infile = file.path(folder, "snps_qc"),
      outfile = outfile
    )
    result <- read.delim(outfile, stringsAsFactors = FALSE)
    stopifnot(all(c("SNP", "CHR", "POS", "PVAL") %in% names(result)),
              nrow(result) > 0L)
    cat("Converged; SNPs:", nrow(result),
        "; GRM variance:", unname(fit$theta[2]), "\n")

    summaries[[paste(name, model_name, sep = "_")]] <- data.frame(
      analysis = name, model = model_name, trait = trait,
      family = likelihood, n_frogs = nrow(pheno),
      n_snps = nrow(result), converged = fit$converged,
      grm_variance = unname(fit$theta[2]),
      stringsAsFactors = FALSE
    )
    rm(fit, result)
    gc(verbose = FALSE)
  }
}

summary <- do.call(rbind, summaries)
write.table(summary, file.path(root, "null_model_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nALL SELECTED GWAS MODELS COMPLETED\n")
print(summary)
