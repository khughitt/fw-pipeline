#!/bin/env Rscript
#
# combines feature weights from multiple sources
#
suppressMessages(library(tidyverse))

if (snakemake@config$dev_mode$enabled) {
    if (!dir.exists('/rda/nih/fw/src/combine')) {
        dir.create('/rda/nih/fw/src/combine', mode = '0755', recursive = TRUE)
    }
  save.image('/rda/nih/fw/src/combine/combine_weights.rda')
}

# feature column name
feat_level <- snakemake@wildcards$feature_level

if (feat_level == 'genes') {
    feat_keys <- 'symbol'
} else if (feat_level == 'gene_sets') {
    feat_keys <- c('collection', 'gene_set')
}

# load individual dataset weights
infiles <- unlist(snakemake@input)
wts_list <- lapply(infiles, read_tsv, col_types = cols())

# TEMP: manually load multiple myeloma-based gene weights for now;
if (length(feat_keys) == 1) {
  # genes
  mm_wts <- read_csv('/data/nih/archive/gene-weights/mm/2019-10-24/combined/gene_weights.csv',
                    col_types = cols()) %>%
    rename(mm_score = score)
} else {
  # gene sets
  mm_wts <- read_csv('/data/nih/archive/gene-weights/mm/2019-10-24/combined/pathway_weights.csv',
                    col_types = cols()) %>%
    rename(mm_score = score)

  # combine <collection, gene_set> cols
  # mm_wts$gene_set <- sprintf("%s_%s", mm_wts$collection, mm_wts$gene_set)

  # mm_wts <- mm_wts %>%
  #   select(-collection)

  # feat_keys <- 'gene_set'
}

# drop mm weights measured in less than 25% of the datasets used for its construction
# and remove uneeded columns
min_meas <- max(mm_wts$num_measurements) * 0.25

mm_wts <- mm_wts %>%
  filter(num_measurements >= min_meas) %>%
  select(feat_keys, mm_score)

wts_list[[length(wts_list) + 1]] <- mm_wts

# reorder weights list so that the entry with the most rows appears first
wts_list <- wts_list[order(-rank(unlist(lapply(wts_list, nrow)), ties.method = 'first'))]

# combine individual weights into a single tibble
wts <- wts_list %>%
  reduce(left_join, by = feat_keys)

# weight_suffix <- sprintf("%s_%s", snakemake@wildcards$collapse_func,
#                          snakemake@wildcards$cor_method)

wts <- wts %>%
  select(feat_keys, mm_score, starts_with('score'))

# scale weights so that each source contributes equally
feat_id_ind <- which(colnames(wts) %in% feat_keys)

unscaled_totals <- colSums(wts[, -feat_id_ind], na.rm = TRUE)

wts[, -feat_id_ind] <- sweep(wts[, -feat_id_ind], 2, unscaled_totals, '/') * 1e6

wts_mat <- wts[, -feat_id_ind]

wts$max_score  <- apply(wts_mat, 1, max, na.rm = TRUE)
wts$mean_score <- apply(wts_mat, 1, mean, na.rm = TRUE)

# combined score
wts <- wts %>%
  arrange(desc(mean_score))

# drop individual scores
wts <- wts %>%
  select(feat_keys, max_score, mean_score)

# store environment
if (snakemake@config$dev_mode$enabled) {
  save.image('/rda/nih/fw/src/combine/combine_weights.rda')
}

write_tsv(wts, snakemake@output[[1]])
