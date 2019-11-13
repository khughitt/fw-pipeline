#
# Measure function enrichment for generate feature weights
#
# TODO: handle duplicate ids in DSigDB and tftargets gene sets
#
suppressMessages(library(GSEABase))
suppressMessages(library(fgsea))
suppressMessages(library(tidyverse))

config <- snakemake@config$gene_sets

# load feature weights
feat_weights <- read_tsv(snakemake@input[[1]], col_types = cols())

# load gene sets
infiles <- file.path(config$include, list.files(config$include))

gene_sets <- lapply(infiles, function(x) {
  res <- geneIds(getGmt(x))
  lapply(res, function(gset) { gset[gset != ''] })
})
names(gene_sets) <- tools::file_path_sans_ext(basename(infiles))

# remove gene set :length suffixes, if present
names(gene_sets) <- sub(':\\d+$', '', names(gene_sets))

# remove any leading or trailing whitespace in gene set names;
# encountered at least one instance of this ("HEK 293 T-rex" in NCI-60 gene set
# collection)
for (gset in names(gene_sets)) {
  names(gene_sets[[gset]]) <- trimws(names(gene_sets[[gset]]))
}

# remove gene set weights, if present
# e.g. "ANXA1,1.0" -> "ANXA1"
for (gset in names(gene_sets)) {
  gene_sets[[gset]] <- lapply(gene_sets[[gset]], function(x) { sub(',\\d+\\.\\d+$', '', x) })
}

# exclude any gene sets which are either too large or too small
for (gset in names(gene_sets)) {
	set_sizes <- lapply(gene_sets[[gset]], length)
	mask <- (set_sizes >= config$min_size) & (set_sizes <= config$max_size)
	gene_sets[[gset]] <- gene_sets[[gset]][mask]
}

fgsea_results <- NULL

# iterate over weight columns
for (i in 2:ncol(feat_weights)) {
  # iterate over gene sets and measure enrichment
  for (gene_set in names(gene_sets)) {
    set.seed(1)

    dat <- deframe(feat_weights[, c(1, i)])

    res <- fgsea(gene_sets[[gene_set]], stats = dat, nperm = config$fgsea_nperm, nproc = 24) %>%
      select(-leadingEdge) %>%
      arrange(pval)

    weight_col <- colnames(feat_weights)[i]
    fgsea_results <- rbind(fgsea_results, cbind(field = weight_col, gene_set, res))
  }
}

write_tsv(fgsea_results, snakemake@output[[1]])

