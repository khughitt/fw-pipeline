
  # match drug and rna sample ids
  sample_ids <- colnames(drug_dat_imputed)
  sample_ids <- sample_ids[sample_ids %in% colnames(rna_dat)]

  drug_dat_imputed <- drug_dat_imputed[, colnames(drug_dat_imputed) %in% sample_ids]

  # create a subsetted version of the rna data with matching samples
  rna_dat_subset <- rna_dat[, colnames(rna_dat) %in% sample_ids]

  # drop any drugs with no variance
  drug_dat_imputed <- drug_dat_imputed[apply(drug_dat_imputed, 1, var) != 0, ]

  # transpose drug data
  drug_dat_imputed <- t(drug_dat_imputed)

  # DEV
  save.image('~/tmp.rda')

  # iterate over correlation methods and measure rna-drug correlations for each
  for (cor_method in cor_config$funcs) {
    drug_gene_cors <- apply(rna_dat_subset, 1, function(gene_expr) {
      drug_cors <- as.numeric(cor(gene_expr, drug_dat_imputed, method = cor_method))
      do.call(cor_config$collapse, list(drug_cors))
    })

    res <- tibble(
      symbol = names(drug_gene_cors),
      cor    = as.numeric(drug_gene_cors)
    )
    colnames(res) <- c('symbol', sprintf("%s_%s_%s", sens_meas, cor_config$collapse, cor_method))

    drug_gene_weights <- drug_gene_weights %>%
      inner_join(res, by = 'symbol')
  }
