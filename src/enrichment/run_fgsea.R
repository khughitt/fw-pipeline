#
# Measure function enrichment for generate feature weights
#
suppressMessages(library(GSEABase))
suppressMessages(library(fgsea))
suppressMessages(library(tidyverse))

config <- snakemake@config$gene_sets

# load feature weights
feat_weights <- read_tsv(snakemake@input[[1]], col_types = cols())

# load gene sets
infiles <- file.path(config$include, list.files(config$include))

gene_set_collections <- lapply(infiles, function(x) {
  res <- geneIds(getGmt(x))
  lapply(res, function(gene_set) {
    gene_set[gene_set != '']
  })
})
names(gene_set_collections) <- tools::file_path_sans_ext(basename(infiles))

# remove gene set :length suffixes, if present
names(gene_set_collections) <- sub(':\\d+$', '', names(gene_set_collections))

# remove any leading or trailing whitespace in gene set names;
# encountered at least one instance of this ("HEK 293 T-rex" in NCI-60 gene set
# collection)
for (collection in names(gene_set_collections)) {
  names(gene_set_collections[[collection]]) <- trimws(names(gene_set_collections[[collection]]))
}

# remove gene set weights, if present
# e.g. "ANXA1,1.0" -> "ANXA1"
for (collection in names(gene_set_collections)) {
  gene_set_collections[[collection]] <- lapply(gene_set_collections[[collection]], function(x) {
    sub(',\\d+\\.\\d+$', '', x)
  })
}

save.image('/rda/nih/fw/run_fgsea.rda')

# exclude any gene sets which are either too large or too small
for (collection in names(gene_set_collections)) {
  set_sizes <- lapply(gene_set_collections[[collection]], length)
  mask <- (set_sizes >= config$min_size) & (set_sizes <= config$max_size)
  gene_set_collections[[collection]] <- gene_set_collections[[collection]][mask]
}

# nov 18, 2019
# work-around for possible corrupted entry in DSigDB GMT (line 6027);
# issue has been reported and manually corrected in the 1.0 currently in use, but
# leaving this code in for now in case it accidentally gets reverted.
ind <- which(names(gene_set_collections[['DSigDB_All']]) == '')

if (length(ind) > 0) {
  gene_set_collections[['DSigDB_All']] <- gene_set_collections[['DSigDB_All']][-ind]
}

fgsea_results <- NULL

# iterate over weight columns
for (i in 2:ncol(feat_weights)) {
  dat <- deframe(feat_weights[, c(1, i)])

  # drop missing values
  dat <- dat[!is.na(dat)]

  weight_col <- colnames(feat_weights)[i]

  # iterate over gene sets and measure enrichment
  for (gene_set in names(gene_set_collections)) {
    set.seed(1)

    res <- fgsea(gene_set_collections[[gene_set]], stats = dat, nperm = config$fgsea_nperm, nproc = 24) %>%
      select(-leadingEdge) %>%
      arrange(pval)

    fgsea_results <- rbind(fgsea_results, cbind(field = weight_col, gene_set, res))
  }
}

# drop entries for which no p-value could be computed
fgsea_results <- fgsea_results[complete.cases(fgsea_results), ]

write_tsv(fgsea_results, snakemake@output[[1]])
