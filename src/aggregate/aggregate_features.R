#!/bin/env Rscript
#
# Aggregates features along known pathways / gene sets
#
suppressMessages(library(GSEABase))
suppressMessages(library(tidyverse))

# load feature data
feat_dat <- read_tsv(snakemake@input[[1]], col_types = cols())

# load gene sets
infile <- file.path(snakemake@config$gene_sets$include, sprintf("%s.gmt", snakemake@wildcards$gene_set))

gene_sets <- geneIds(getGmt(infile))
gene_sets <- gene_sets[gene_sets != '']

# remove any leading or trailing whitespace in gene set names;
# encountered at least one instance of this ("HEK 293 T-rex" in NCI-60 gene set
# collection)
names(gene_sets) <- trimws(names(gene_sets))

# remove gene set weights, if present
# e.g. "ANXA1,1.0" -> "ANXA1"
gene_sets <- lapply(gene_sets, function(x) {
  sub(',\\d+\\.\\d+$', '', x)
})

# exclude any gene sets which are either too large or too small
set_sizes <- lapply(gene_sets, length)
mask <- (set_sizes >= snakemake@config$gene_sets$min_size) & (set_sizes <= snakemake@config$gene_sets$max_size)
gene_sets <- gene_sets[mask]

# create a pathway / gene set-projected version of feature data; collection and gene set
# information are stored in a separate data.frame to preserve numeric name of feature
# data
aggregated_feat_id_cols <- NULL
aggregated_feat_dat <- NULL

for (gene_set in names(gene_sets)) {
  feat_dat_subset <- feat_dat[feat_dat$symbol %in% gene_sets[[gene_set]], ]

  # skip any genes sets with too few genes present in data
  if (nrow(feat_dat_subset) < snakemake@config$gene_sets$min_size) {
    next
  }

  dat <- apply(feat_dat_subset[, -1], 2, snakemake@wildcards$agg_func, na.rm = TRUE)

  aggregated_feat_id_cols <- rbind(aggregated_feat_id_cols, c(snakemake@wildcards$gene_set, gene_set))
  aggregated_feat_dat <- rbind(aggregated_feat_dat, dat)
}

# combine aggregated data and id columns
aggregated_feat_dat <- cbind(as.data.frame(aggregated_feat_id_cols),
                             aggregated_feat_dat)

colnames(aggregated_feat_dat)[1:2] <- c('collection', 'gene_set')

# drop any rows with zero variance (uninformative)
mask <- apply(aggregated_feat_dat[, 3:ncol(aggregated_feat_dat)], 1, var, na.rm = TRUE) > 0
aggregated_feat_dat <- aggregated_feat_dat[mask, ]

write_tsv(aggregated_feat_dat, snakemake@output[[1]])
