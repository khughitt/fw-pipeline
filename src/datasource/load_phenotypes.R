#
# snakemake phenotype parsing script
#
# chooses script to execute baed on the datasource type.
#
options(stringsAsFactors = FALSE)
set.seed(1)

config <- snakemake@params$config

# detect phenotype and data source type
phenotype   <- snakemake@wildcards$phenotype

print(phenotype)
print(config$phenotype)

pheno_config <- config$phenotypes[[phenotype]]
data_source <- pheno_config$data_source

# pharmacogx
if (data_source == 'pharmacogx') {
  source('src/datasource/pharmacogx/load_phenotypes.R')
  load_pharmacogx_phenotype(config, snakemake@output[[1]])
}
