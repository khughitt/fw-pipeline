#/usr/bin Rscript
#
# Generate clean version of PharmacoGx drug data
#
suppressMessages(library(annotables))
suppressMessages(library(Biobase))
suppressMessages(library(PharmacoGx))
suppressMessages(library(VIM))
suppressMessages(library(tidyverse))

load_pharmacogx_phenotype <- function(config, output_file) {
  pset_id <- config$pset

  supported_psets <- c('GDSC1000', 'CCLE_2013')

  if (!pset_id %in% supported_psets) {
    stop(sprintf("Unsupported Pharmacoset specified: %s!", pset_id))
  }

  pset_rda <- file.path(config$cache_dir, paste0(pset_id, '.RData'))

  if (!file.exists(pset_rda)) {
    pset <- downloadPSet(pset_id, saveDir = config$cache_dir)
  } else {
    pset <- get(load(pset_rda))
  }

  # iterate over phenotype datasets and generate "clean" versions of each
  drug_dat <- summarizeSensitivityProfiles(pset, phenotype, 
                                          summary.stat = config$drug[[phenotype]]$summary_stat)
  drug_config <- config$phenotypes[[phenotype]]

  # clip extreme values
  if ('clip' %in% names(drug_config)) {
    clip_lower <- as.numeric(drug_config$clip$min_val)
    clip_upper <- as.numeric(drug_config$clip$max_val)

    drug_dat[drug_dat < clip_lower] <- clip_lower
    drug_dat[drug_dat > clip_upper] <- clip_upper
  }

  # log-transform sensitivity scores
  if (drug_config$log_transform) {
    drug_dat <- log1p(drug_dat)
  }

  # remove samples with too many missing values
  if ('filter' %in% names(drug_config)) {
    row_num_na <- apply(drug_dat, 1, function(x) { sum(is.na(x)) })
    col_num_na <- apply(drug_dat, 2, function(x) { sum(is.na(x)) })

    drug_data_col_max_na <- round((nrow(drug_dat) - 1) * drug_config$filter$col_max_na)
    drug_data_row_max_na <- round((ncol(drug_dat) - 1) * drug_config$filter$row_max_na)

    # filter samples / drugs with too many missing values
    drug_dat <- drug_dat[, col_num_na <= drug_data_col_max_na]
    drug_dat <- drug_dat[row_num_na <= drug_data_row_max_na, ]
  }

  # impute remaining missing values
  drug_dat_imputed <- as.matrix(kNN(t(drug_dat), k = drug_config$impute$k)[, 1:nrow(drug_dat)])
  rownames(drug_dat_imputed) <- colnames(drug_dat)

  # drop any drugs with no variance
  drug_dat_imputed <- drug_dat_imputed[, apply(drug_dat_imputed, 2, var) != 0]

  # transpose back to original orientation and store
  drug_dat_imputed <- t(drug_dat_imputed)

  drug_dat_imputed %>%
    as.data.frame() %>%
    rownames_to_column('drug') %>%
    write_tsv(output_file)
}
