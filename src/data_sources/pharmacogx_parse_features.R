#/usr/bin Rscript
#
# Generate clean version of PharmacoGx feature data
#
suppressMessages(library(annotables))
suppressMessages(library(Biobase))
suppressMessages(library(PharmacoGx))
suppressMessages(library(tidyverse))

options(stringsAsFactors = FALSE)
set.seed(1)

#
# Setup
#
config <- snakemake@params$config
pset_id <- config$pset

#cor_config <- snakemake@config$methods$correlation

supported_psets <- c('GDSC1000')

if (!pset_id %in% supported_psets) {
  stop(sprintf("Unsupported Pharmacoset specified: %s!", pset_id))
}

pset_rda <- file.path(config$rda_dir, paste0(pset_id, '.RData'))

if (!file.exists(pset_rda)) {
  pset <- downloadPSet(pset_id, saveDir = config$rda_dir)
} else {
  pset <- get(load(pset_rda))
}

# iterate over feature datasets and generate "clean" versions of each
for (mdatatype in names(config$features)) {
  eset <- summarizeMolecularProfiles(pset, mDataType = mdatatype)
  dat <- bind_cols(symbol = fData(eset)$Symbol, as.data.frame(exprs(eset)))

  # retrieve missing gene symbols from annotables
  missing_ind <- which(is.na(dat$symbol))
  missing_ensgene <- fData(eset)$Ensembl[missing_ind]

  dat$symbol[missing_ind] <- grch38$symbol[match(missing_ensgene, grch38$ensgene)]

  # drop rna entries that could not be mapped to gene symbols
  dat <- dat[!is.na(dat$symbol), ]

  # drop any samples that are completely missing gene expression data
  all_missing <- function(x) {
    sum(is.na(x)) == length(x)
  }
  dat <- dat[, !apply(dat, 2, all_missing)]

  # collapse multi-mapped entries and store result
  dat <- dat %>%
    group_by(symbol) %>%
    summarize_all(mean) %>%
    ungroup %>%
    write_tsv(snakemake@output[[1]])
}
