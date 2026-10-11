library(GMMAT)


pheno <- read.delim(
  file.path(gwas_dir, "gwas_input.tsv")
)

K <- as.matrix(
  read.csv(
    file.path(gwas_dir, "grm_892.csv"),
    header = FALSE
  )
)

grm_ids <- read.delim(
  file.path(gwas_dir, "grm_892_ids.tsv")
)

stopifnot(identical(pheno$IID, grm_ids$IID))

rownames(K) <- grm_ids$IID
colnames(K) <- grm_ids$IID

# Initial model: GRM covariance, intercept only.
model_grm <- glmmkin(
  infected_2plus ~ 1,
  data = pheno,
  kins = K,
  id = "IID",
  family = binomial(link = "logit")
)