#!/bin/env Rscript
#
# combines feature weights from multiple sources
#
suppressMessages(library(tidyverse))

# load individual dataset weights
infiles <- unlist(snakemake@input)
wts_list <- lapply(infiles, read_tsv, col_types = cols())

# add prefix to dataset columns;
# disabled for now -- too messy and not needed 
#str_split(sub(snakemake@config$output_dir, '', infiles), '/')

# TEMP: manually load multiple myeloma-based gene weights for now;
mm_wts <- read_csv('/data/nih/gene-weights/mm/2019-10-24/combined/gene_weights.csv') %>%
  rename(mm_score = score)

# drop mm weights measured in less than 25% of the datasets used for its construction
# and remove uneeded columns
min_meas <- max(mm_wts$num_measurements) * 0.25

mm_wts <- mm_wts %>%
  filter(num_measurements >= min_meas) %>%
  select(symbol, mm_score)

wts_list[[length(wts_list) + 1]] <- mm_wts

# reorder list so that entry with the most rows appears first
wts_list <- wts_list[rank(-unlist(lapply(wts_list, nrow)))]

# combine individual weights into a single tibble
wts <- wts_list %>%
  reduce(left_join, by = "symbol")

# temp: hard-coded list of columns of interest (mm/ccle/gdsc)
wts <- wts %>%
  select(symbol, mm_score, ends_with('max_abs_spearman'))

# replace missing values with column medians
mask <- sapply(wts, function(x) { any(is.na(x)) })

for (field in colnames(wts)[mask]) {
  vals <- pull(wts, field)
  vals[is.na(vals)] <- median(vals, na.rm = TRUE)
  wts[, field] <- vals
}

# convert scores to ranks
wts <- wts %>%
  mutate_if(is.numeric, function(x) { rank(-x) })

wts$mean_rank <- apply(wts[, -1], 1, function(x) { mean(x, na.rm = TRUE) })
wts$min_rank  <- apply(wts[, -1], 1, min, na.rm = TRUE )

# combined score: min(min_rank, rank(mean_rank))
wts <- wts %>%
  mutate(rank_mean_rank = rank(mean_rank)) %>%
  mutate(combined_rank = pmin(min_rank, rank_mean_rank)) %>%
  select(symbol, min_rank, mean_rank, combined_rank) %>%
  arrange(combined_rank)

write_tsv(wts, snakemake@output[[1]])
