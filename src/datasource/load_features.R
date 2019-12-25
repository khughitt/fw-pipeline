#
# snakemake feature parsing script
#
suppressMessages(library(tidyverse))

options(stringsAsFactors = FALSE)
set.seed(1)

# data source config
config <- snakemake@params$config

# detect feature and data source type
feature_type <- snakemake@wildcards$feature_type
feature_level <- snakemake@wildcards$feature_level

# feature config
feature_config <- config$features[[feature_type]][[feature_level]]

data_source <- feature_config$data_source

save.image('~/tmp-parse-features.rda')

# pharmacogx
if (data_source == 'pharmacogx') {
  source('src/datasource/pharmacogx/load_features.R')
  load_pharmacogx_feature(config, feature_config$mdatatype, snakemake@output[[1]])
} else if (data_source == 'tsv') {
  file.symlink(feature_config$path, snakemake@output[[1]])
}
