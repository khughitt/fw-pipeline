#
# Generate clean version of PharmacoGx feature data
#
suppressMessages(library(annotables))
suppressMessages(library(Biobase))
suppressMessages(library(PharmacoGx))

#
# Setup
#
load_pharmacogx_feature <- function(config, mdatatype, output_file) {
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

  eset <- summarizeMolecularProfiles(pset, mDataType = mdatatype)

  # convert expressionset to a tibble
  if ('symbol' %in% colnames(fData(eset))) {
    gene_symbols <- fData(eset)$symbol
  } else if ('Symbol' %in% colnames(fData(eset))) {
    gene_symbols <- fData(eset)$Symbol
  }
  dat <- bind_cols(symbol = gene_symbols, as.data.frame(exprs(eset)))

  # retrieve missing gene symbols from annotables
  if ('EnsemblGeneId' %in% colnames(fData(eset))) {
    # GDSC, GDSC1000 (GRCh37 / ensgene)
    missing_ind <- which(is.na(dat$symbol))
    missing_ensgene <- fData(eset)$EnsemblGeneId[missing_ind]

    dat$symbol[missing_ind] <- grch38$symbol[match(missing_ensgene, grch38$ensgene)]
  } else if ('GENEID' %in% colnames(fData(eset))) {
    # CCLE (GRCh38 / entrez)
    missing_ind <- which(is.na(dat$symbol))
    missing_entrez <- fData(eset)$GENEID[missing_ind]

    dat$symbol[missing_ind] <- grch37$symbol[match(missing_entrez, grch37$entrez)]
  }

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
    write_tsv(output_file)
}
