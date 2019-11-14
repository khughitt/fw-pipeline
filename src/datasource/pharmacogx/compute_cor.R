#
# Compute PharmacoGx drug-feature correlations
#
suppressMessages(library(tidyverse))

# get correlation settings
config <- snakemake@config$methods$correlation

# load processed PharmacoGx drug and molecular feature data
feat_dat <- read_tsv(snakemake@input$features)
drug_dat <- read_tsv(snakemake@input$phenotype)

# match drug and feature sample ids
sample_ids <- colnames(drug_dat)[-1]
sample_ids <- sort(sample_ids[sample_ids %in% colnames(feat_dat)[-1]])

feat_id_col <- colnames(feat_dat)[1]
drug_id_col <- colnames(drug_dat)[1]

feat_dat <- feat_dat[, c(feat_id_col, sample_ids)]
drug_dat <- drug_dat[, c(drug_id_col, sample_ids)]

# transpose drug data to avoid having to do multiple times later on
drug_dat <- t(drug_dat[, -1])

# iterate over correlation methods and measure feature-drug correlations for each
drug_feat_weights <- NULL

# helper functions; will be moved to external file..
max_abs <- function(x) {
  max(abs(x), na.rm = TRUE)
}

for (cor_method in config$funcs) {
  for (collapse_func in config$collapse) {
    drug_feat_cors <- apply(feat_dat[, -1], 1, function(feat_vals) {
      drug_cors <- as.numeric(cor(feat_vals, drug_dat[, -1], method = cor_method))
      do.call(collapse_func, list(drug_cors))
    })

    res <- tibble(
      symbol = pull(feat_dat, 1),
      cor    = as.numeric(drug_feat_cors)
    )
    colnames(res) <- c(feat_id_col, sprintf("%s_%s_%s", snakemake@wildcards$phenotype, collapse_func, cor_method))

    if (is.null(drug_feat_weights)) {
      drug_feat_weights <- res
    } else {
      drug_feat_weights <- drug_feat_weights %>%
        inner_join(res, by = feat_id_col)
    }
  }
}

# save result
write_tsv(drug_feat_weights, snakemake@output[[1]])

