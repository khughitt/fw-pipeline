#!/bin/env Rscript
#
# combines feature weights from multiple sources
#
suppressMessages(library(tidyverse))

# load individual dataset weights
infiles <- unlist(snakemake@input)
wts_list <- lapply(infiles, read_tsv, col_types = cols())

# TEMP: manually load multiple myeloma-based gene weights for now;
mm_wts <- read_csv('/data/nih/gene-weights/mm/2019-10-24/combined/gene_weights.csv',
                   col_types = cols()) %>%
  rename(mm_score = score)

# drop mm weights measured in less than 25% of the datasets used for its construction
# and remove uneeded columns
min_meas <- max(mm_wts$num_measurements) * 0.25

mm_wts <- mm_wts %>%
  filter(num_measurements >= min_meas) %>%
  select(symbol, mm_score)

wts_list[[length(wts_list) + 1]] <- mm_wts

# reorder weights list so that the entry with the most rows appears first
wts_list <- wts_list[order(-rank(unlist(lapply(wts_list, nrow)), ties.method = 'first'))]

# combine individual weights into a single tibble
wts <- wts_list %>%
  reduce(left_join, by = "symbol")

weight_suffix <- sprintf("%s_%s", snakemake@wildcards$collapse_func,
                         snakemake@wildcards$cor_method)

wts <- wts %>%
  select(symbol, mm_score, ends_with(weight_suffix))

# scale weights so that each source contributes equally
wts[, -1] <- scale(wts[, -1])

wts_mat <- wts[, -1]

wts$max_score  <- apply(wts_mat, 1, max, na.rm = TRUE)
wts$mean_score <- apply(wts_mat, 1, mean, na.rm = TRUE)

# combined score
wts <- wts %>%
  arrange(desc(mean_score))

# drop individual scores
wts <- wts %>%
  select(symbol, max_score, mean_score)

# store environment
if (snakemake@config$dev_mode$enabled) {
  save.image('/rda/nih/fw/combine_weights.rda')
}

write_tsv(wts, snakemake@output[[1]])
