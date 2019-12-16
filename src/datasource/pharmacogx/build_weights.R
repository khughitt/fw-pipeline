#
# Compute PharmacoGx drug-feature correlations
#
suppressMessages(library(tidyverse))

# get correlation settings
correlation_params <- snakemake@config$combine$correlation

# load processed PharmacoGx drug and molecular feature data
feat_dat <- read_tsv(snakemake@input$features)
drug_dat <- read_tsv(snakemake@input$phenotype)

# determine feature sample indices
if (snakemake@params$feature_level == 'genes') {
  # sample_id
  feat_id_ind <- 1
} else {
  # collection, gene_set
  feat_id_ind <- 1:2
}

# match drug and feature sample ids
sample_ids <- colnames(drug_dat)[-1]
sample_ids <- sort(sample_ids[sample_ids %in% colnames(feat_dat)[-feat_id_ind]])

feat_id_cols <- colnames(feat_dat)[feat_id_ind]
drug_id_col <- colnames(drug_dat)[1]

feat_dat <- feat_dat[, c(feat_id_cols, sample_ids)]
drug_dat <- drug_dat[, c(drug_id_col, sample_ids)]

feat_id_ind <- colnames(feat_dat) %in% feat_id_cols
feat_ind    <- !feat_id_ind

# transpose drug data now to avoid having to do multiple times later on
drug_dat <- t(drug_dat[, -1])

# iterate over correlation methods and measure feature-drug correlations for each
drug_feat_weights <- NULL

if (snakemake@config$dev_mode$enabled) {
  save.image(file.path(snakemake@config$dev_mode$rda_dir, 'build_weights.rda'))
}

# helper functions; will be moved to external file..
maxabs <- function(x) {
  max(abs(x), na.rm = TRUE)
}

for (cor_method in correlation_params$methods) {
  for (collapse_func in correlation_params$collapse_funcs) {
    drug_feat_cors <- apply(feat_dat[, feat_ind], 1, function(feat_vals) {
      drug_cors <- as.numeric(cor(feat_vals, drug_dat, method = cor_method))
      do.call(collapse_func, list(drug_cors))
    })

    if (snakemake@params$feature_level == 'genes') {
      # gene result
      res <- tibble(
        symbol = pull(feat_dat, 1),
        cor    = as.numeric(drug_feat_cors)
      )
    } else {
      # gene set result
      res <- tibble(
        collection = pull(feat_dat, 1),
        gene_set   = pull(feat_dat, 2),
        cor        = as.numeric(drug_feat_cors)
      )
    }

    cname <- sprintf("%s_%s_%s_%s", snakemake@wildcards$dataset,
                     snakemake@wildcards$phenotype, collapse_func, cor_method)
    colnames(res) <- c(feat_id_cols, cname)

    if (is.null(drug_feat_weights)) {
      drug_feat_weights <- res
    } else {
      drug_feat_weights <- drug_feat_weights %>%
        inner_join(res, by = feat_id_cols)
    }
  }
}

if (snakemake@config$dev_mode$enabled) {
  save.image(file.path(snakemake@config$dev_mode$rda_dir, 'build_weights.rda'))
}

# save result
write_tsv(drug_feat_weights, snakemake@output[[1]])
