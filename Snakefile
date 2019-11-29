"""
Feature Weights Pipeline
V. Keith Hughitt
"""
import glob
import os
import pandas as pd
import yaml
from os.path import join
from pathlib import Path

# load data source configs
configs = {}

for cfg in glob.glob('data_sources/*.yml'):
    with open(cfg, 'r') as fp:
        source_config = yaml.safe_load(fp)
        configs[source_config['name']] = source_config

# get a list of gene sets to be analyzed
gene_set_names = [x.replace('.gmt', '') for x in os.listdir(config['gene_sets']['include'])]

# Generate lists of output filepath components
datasets   = []
versions   = []
features   = []
phenotypes = []

for dataset_id in configs:
    for feature_type in configs[dataset_id]['features']:
        for phenotype in configs[dataset_id]['phenotypes']:
            datasets.append(dataset_id)
            versions.append(configs[dataset_id]['version'])
            features.append(feature_type)
            phenotypes.append(phenotype)

# same as above, but for gene set outputs
gene_set_datasets   = []
gene_set_versions   = []
gene_set_features   = []
gene_set_phenotypes = []

for dataset_id in configs:
    for feature_type in configs[dataset_id]['features']:
        for phenotype in configs[dataset_id]['phenotypes']:
            gene_set_datasets.append(dataset_id)
            gene_set_versions.append(configs[dataset_id]['version'])
            gene_set_features.append(feature_type)
            gene_set_phenotypes.append(phenotype)

# weight combinations
cor_methods = config['combine']['correlation']['methods']
collapse_funcs = config['combine']['correlation']['collapse_funcs']

# create directory to store environment states (debug-mode)
if config['dev_mode']['enabled']:
    os.makedirs(config['dev_mode']['rda_dir'], mode = 755, exist_ok=True)

rule summarize_combined_weights:
    input: 
        indiv_weights = expand(join(config['output_dir'], 'weights/individual/{dataset}/{version}/genes/{feature}/{phenotype}.tsv.gz'), zip, dataset=datasets, version=versions, feature=features, phenotype=phenotypes),
        indiv_fgsea = expand(join(config['output_dir'], 'fgsea/individual/{dataset}/{version}/genes/{feature}/{phenotype}.tsv.gz'), zip, dataset=datasets, version=versions, feature=features, phenotype=phenotypes),
        combined_weights = expand(join(config['output_dir'], 'weights/combined/genes/combined_weights_{cor_method}_{collapse_func}.tsv.gz'), cor_method=cor_methods, collapse_func=collapse_funcs),
        combined_fgsea = expand(join(config['output_dir'], 'fgsea/combined/genes/combined_weights_{cor_method}_{collapse_func}_fgsea.tsv.gz'), cor_method=cor_methods, collapse_func=collapse_funcs)
    output:
        join(config['report_dir'], config['version'], 'combined_weights_genes_summary.html')
    script:
        'reports/combined_weights_summary.Rmd'

rule summarize_combined_aggregated_weights:
    input: 
        indiv_weights = expand(join(config['output_dir'], 'weights/individual/{dataset}/{version}/gene_sets/{feature}/{phenotype}.tsv.gz'), zip, dataset=datasets, version=versions, feature=features, phenotype=phenotypes),
        indiv_fgsea = expand(join(config['output_dir'], 'fgsea/individual/{dataset}/{version}/gene_sets/{feature}/{phenotype}.tsv.gz'), zip, dataset=datasets, version=versions, feature=features, phenotype=phenotypes),
        combined_weights = expand(join(config['output_dir'], 'weights/combined/gene_sets/combined_weights_{cor_method}_{collapse_func}.tsv.gz'), cor_method=cor_methods, collapse_func=collapse_funcs),
        combined_fgsea = expand(join(config['output_dir'], 'fgsea/combined/gene_sets/combined_weights_{cor_method}_{collapse_func}_fgsea.tsv.gz'), cor_method=cor_methods, collapse_func=collapse_funcs)
    output:
        join(config['report_dir'], config['version'], 'combined_weights_gene_sets_summary.html')
    script:
        'reports/gene_sets/combined_weights_summary.Rmd'

rule run_fgsea_combined_individual_weights:
    input: 
        join(config['output_dir'], 'weights/combined/genes/combined_weights_{cor_method}_{collapse_func}.tsv.gz')
    output:
        join(config['output_dir'], 'fgsea/combined/genes/combined_weights_{cor_method}_{collapse_func}_fgsea.tsv.gz')
    script: 'src/enrichment/run_fgsea.R'

rule run_fgsea_combined_aggregated_weights:
    input: 
        #expand(join(config['output_dir'], 'weights/combined_weights_{cor_method}_{collapse_func}s.tsv.gz'), 
        #       cor_method=cor_methods, collapse_func=collapse_funcs)
        join(config['output_dir'], 'weights/combined/gene_sets/combined_weights_{cor_method}_{collapse_func}.tsv.gz')
    output:
        join(config['output_dir'], 'fgsea/combined/gene_sets/combined_weights_{cor_method}_{collapse_func}_fgsea.tsv.gz')
    script: 'src/enrichment/run_fgsea.R'

rule run_fgsea_individual_weights:
    input: 
        join(config['output_dir'], 'weights/individual/{dataset}/{version}/genes/{feature}/{phenotype}.tsv.gz') 
    output:
        join(config['output_dir'], 'fgsea/individual/{dataset}/{version}/genes/{feature}/{phenotype}.tsv.gz')
    script: 'src/enrichment/run_fgsea.R'

rule run_fgsea_aggregated_weights:
    input: 
        join(config['output_dir'], 'weights/individual/{dataset}/{version}/gene_sets/{feature}/{phenotype}.tsv.gz') 
    output:
        join(config['output_dir'], 'fgsea/individual/{dataset}/{version}/gene_sets/{feature}/{phenotype}.tsv.gz')
    script: 'src/enrichment/run_fgsea.R'

rule combine_individual_weights:
    input: 
        expand(join(config['output_dir'], 'weights/individual/{dataset}/{version}/genes/{feature}/{phenotype}.tsv.gz'),
               zip,
               dataset=datasets, version=versions, feature=features, phenotype=phenotypes)
    output:
        join(config['output_dir'], 'weights/combined/genes/combined_weights_{cor_method}_{collapse_func}.tsv.gz')
    script: 'src/combine/combine_weights.R'

rule combine_aggregated_weights:
    input: 
        expand(join(config['output_dir'], 'weights/individual/{dataset}/{version}/gene_sets/{feature}/{phenotype}.tsv.gz'),
               zip,
               #dataset=datasets, version=versions, feature=features, phenotype=phenotypes)
        dataset=gene_set_datasets, version=gene_set_versions, feature=gene_set_features, phenotype=gene_set_phenotypes)
    output:
        join(config['output_dir'], 'weights/combined/gene_sets/combined_weights_{cor_method}_{collapse_func}.tsv.gz')
    script: 'src/combine/combine_weights.R'

rule build_individual_weights:
    input:
        features=join(config['output_dir'], 'datasets/{dataset}/{version}/features/genes/{feature}.tsv.gz'),
        phenotype=join(config['output_dir'], 'datasets/{dataset}/{version}/phenotypes/{phenotype}.tsv.gz')
    output:
        join(config['output_dir'], 'weights/individual/{dataset}/{version}/genes/{feature}/{phenotype}.tsv.gz')
    script: 'src/datasource/pharmacogx/build_weights.R'

rule build_aggregated_weights:
    input:
        features=join(config['output_dir'], 'datasets/{dataset}/{version}/features/gene_sets/{feature}_gene_sets.tsv.gz'),
        phenotype=join(config['output_dir'], 'datasets/{dataset}/{version}/phenotypes/{phenotype}.tsv.gz')
    output:
        join(config['output_dir'], 'weights/individual/{dataset}/{version}/gene_sets/{feature}/{phenotype}.tsv.gz')
    script: 'src/datasource/pharmacogx/build_weights.R'

rule merge_aggregated_features:
    input:
        expand(join(config['output_dir'], 'datasets/{{dataset}}/{{version}}/features/gene_sets/{gene_set}/{{feature}}.tsv.gz'), gene_set=gene_set_names)
    output:
        join(config['output_dir'], 'datasets/{dataset}/{version}/features/gene_sets/{feature}_gene_sets.tsv.gz')

rule aggregate_features:
    input:
        join(config['output_dir'], 'datasets/{dataset}/{version}/features/genes/{feature}.tsv.gz')
    output:
        join(config['output_dir'], 'datasets/{dataset}/{version}/features/gene_sets/{gene_set}/{feature}.tsv.gz')
    script: 'src/aggregate/aggregate_features.R'

rule parse_features:
    output:
        join(config['output_dir'], 'datasets/{dataset}/{version}/features/genes/{feature}.tsv.gz')
    params:
        config=lambda wildcards, output: configs[wildcards.dataset]
    script: 'src/datasource/pharmacogx/parse_features.R'

rule parse_phenotypes:
    output:
        join(config['output_dir'], 'datasets/{dataset}/{version}/phenotypes/{phenotype}.tsv.gz')
    params:
        config=lambda wildcards, output: configs[wildcards.dataset]
    script: 'src/datasource/pharmacogx/parse_phenotypes.R'

