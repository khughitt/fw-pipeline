#
# Compute feature-phenotype correlations
#
suppressMessages(library(tidyverse))

# get datasource settings
data_source_config <- snakemake@config$data_sources[[snakemake@wildcards$data_source]]

# get correlation settings
correlation_params <- snakemake@config$combine$correlation

# load processed feature / phenotype data 
feat_dat <- read_tsv(snakemake@input$features, col_types = cols())
pheno_dat <- read_tsv(snakemake@input$phenotype, col_types = cols())

# determine feature sample indices
feat_type <- snakemake@wildcards$feature_type
feat_level <- snakemake@wildcards$feature_level

feat_keys <- data_source_config$features[[feat_type]][[feat_level]]$keys
feat_keys_ind <- which(colnames(feat_dat) == feat_keys)

# match phenotype and feature sample ids
sample_ids <- colnames(pheno_dat)[-1]
sample_ids <- sort(sample_ids[sample_ids %in% colnames(feat_dat)[-feat_keys_ind]])

pheno_key <- colnames(pheno_dat)[1]

feat_dat <- feat_dat[, c(feat_keys, sample_ids)]
pheno_dat <- pheno_dat[, c(pheno_key, sample_ids)]

# transpose phenotype data now to avoid having to do multiple times later on
pheno_dat <- t(pheno_dat[, -1])

if (snakemake@config$dev_mode$enabled) {
  save.image(file.path(snakemake@config$dev_mode$rda_dir, 'build_weights.rda'))
}

# helper functions; will be moved to external file..
maxabs <- function(x) {
  max(abs(x), na.rm = TRUE)
}

# iterate over correlation methods and measure feature-phenotype correlations for each
pheno_feat_weights <- NULL

for (cor_method in correlation_params$methods) {
  for (collapse_func in correlation_params$collapse_funcs) {
    pheno_feat_cors <- apply(feat_dat[, -feat_keys_ind], 1, function(feat_vals) {
      pheno_cors <- as.numeric(cor(feat_vals, pheno_dat, method = cor_method))
      do.call(collapse_func, list(pheno_cors))
    })

    res <- bind_cols(
      select(feat_dat, feat_keys),
      cor    = as.numeric(pheno_feat_cors)
    )

    cname <- sprintf("%s_%s_%s_%s", snakemake@wildcards$data_source,
                     snakemake@wildcards$phenotype, collapse_func, cor_method)
    colnames(res) <- c(feat_keys, cname)

    if (is.null(pheno_feat_weights)) {
      pheno_feat_weights <- res
    } else {
      pheno_feat_weights <- pheno_feat_weights %>%
        inner_join(res, by = feat_keys)
    }
  }
}

if (snakemake@config$dev_mode$enabled) {
  save.image(file.path(snakemake@config$dev_mode$rda_dir, 'build_weights.rda'))
}

# save result
write_tsv(pheno_feat_weights, snakemake@output[[1]])
