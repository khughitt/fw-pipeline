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

# output directory
output_dir = join(config['output_dir'], config['version'])
report_dir = join(config['report_dir'], config['version'])

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
features   = []
phenotypes = []

for dataset_id in configs:
    for feature_type in configs[dataset_id]['features']:
        for phenotype in configs[dataset_id]['phenotypes']:
            datasets.append(dataset_id)
            features.append(feature_type)
            phenotypes.append(phenotype)

# same as above, but for aggregated features
agg_datasets   = []
agg_features   = []
agg_phenotypes = []
agg_funcs      = []

for dataset_id in configs:
    for feature_type in configs[dataset_id]['features']:
        for agg_func in config['gene_sets']['aggregation_funcs']:
            for phenotype in configs[dataset_id]['phenotypes']:
                agg_datasets.append(dataset_id)
                agg_features.append(feature_type)
                agg_funcs.append(agg_func)
                agg_phenotypes.append(phenotype)

# weight combinations
cor_methods = config['combine']['correlation']['methods']
collapse_funcs = config['combine']['correlation']['collapse_funcs']

# create directory to store environment states (debug-mode)
if config['dev_mode']['enabled']:
    os.makedirs(config['dev_mode']['rda_dir'], mode = 755, exist_ok=True)

rule summarize_combined_aggregated_weights:
    input: 
        indiv_weights = expand(join(output_dir, 'weights/individual/{dataset}/gene_sets/{agg_func}/{feature}/{phenotype}.tsv.gz'), zip, dataset=datasets, agg_func=agg_funcs, feature=features, phenotype=phenotypes),
        indiv_fgsea = expand(join(output_dir, 'fgsea/individual/{dataset}/gene_sets/{agg_func}/{feature}/{phenotype}.tsv.gz'), zip, dataset=datasets, agg_func=agg_funcs, feature=features, phenotype=phenotypes),
        combined_weights = expand(join(output_dir, 'weights/combined/gene_sets/{agg_func}/combined_weights_{cor_method}_{collapse_func}.tsv.gz'), agg_func=agg_funcs, cor_method=cor_methods, collapse_func=collapse_funcs),
        combined_fgsea = expand(join(output_dir, 'fgsea/combined/gene_sets/{agg_func}/combined_weights_{cor_method}_{collapse_func}_fgsea.tsv.gz'), agg_func=agg_funcs, cor_method=cor_methods, collapse_func=collapse_funcs)
    output:
        join(report_dir, 'combined_aggregated_weights_summary.html')
    script:
        'reports/combined_aggregated_weights_summary.Rmd'

rule summarize_combined_weights:
    input: 
        indiv_weights = expand(join(output_dir, 'weights/individual/{dataset}/genes/{feature}/{phenotype}.tsv.gz'), zip, dataset=datasets, feature=features, phenotype=phenotypes),
        indiv_fgsea = expand(join(output_dir, 'fgsea/individual/{dataset}/genes/{feature}/{phenotype}.tsv.gz'), zip, dataset=datasets, feature=features, phenotype=phenotypes),
        combined_weights = expand(join(output_dir, 'weights/combined/genes/combined_weights_{cor_method}_{collapse_func}.tsv.gz'), cor_method=cor_methods, collapse_func=collapse_funcs),
        combined_fgsea = expand(join(output_dir, 'fgsea/combined/genes/combined_weights_{cor_method}_{collapse_func}_fgsea.tsv.gz'), cor_method=cor_methods, collapse_func=collapse_funcs)
    output:
        join(report_dir, 'combined_weights_genes_summary.html')
    script:
        'reports/combined_weights_summary.Rmd'

rule run_fgsea_combined_individual_weights:
    input: 
        join(output_dir, 'weights/combined/genes/combined_weights_{cor_method}_{collapse_func}.tsv.gz')
    output:
        join(output_dir, 'fgsea/combined/genes/combined_weights_{cor_method}_{collapse_func}_fgsea.tsv.gz')
    script: 'src/enrichment/run_fgsea.R'

rule run_fgsea_combined_aggregated_weights:
    input: 
        #expand(join(output_dir, 'weights/combined_weights_{cor_method}_{collapse_func}s.tsv.gz'), 
        #       cor_method=cor_methods, collapse_func=collapse_funcs)
        join(output_dir, 'weights/combined/gene_sets/{agg_func}/combined_weights_{cor_method}_{collapse_func}.tsv.gz')
    output:
        join(output_dir, 'fgsea/combined/gene_sets/{agg_func}/combined_weights_{cor_method}_{collapse_func}_fgsea.tsv.gz')
    script: 'src/enrichment/run_fgsea.R'

rule run_fgsea_individual_weights:
    input: 
        join(output_dir, 'weights/individual/{dataset}/genes/{feature}/{phenotype}.tsv.gz') 
    output:
        join(output_dir, 'fgsea/individual/{dataset}/genes/{feature}/{phenotype}.tsv.gz')
    script: 'src/enrichment/run_fgsea.R'

rule run_fgsea_aggregated_weights:
    input: 
        join(output_dir, 'weights/individual/{dataset}/gene_sets/{agg_func}/{feature}/{phenotype}.tsv.gz') 
    output:
        join(output_dir, 'fgsea/individual/{dataset}/gene_sets/{agg_func}/{feature}/{phenotype}.tsv.gz')
    script: 'src/enrichment/run_fgsea.R'

rule combine_individual_weights:
    input: 
        expand(join(output_dir, 'weights/individual/{dataset}/genes/{feature}/{phenotype}.tsv.gz'),
               zip,
               dataset=datasets, feature=features, phenotype=phenotypes)
    output:
        join(output_dir, 'weights/combined/genes/combined_weights_{cor_method}_{collapse_func}.tsv.gz')
    script: 'src/combine/combine_weights.R'

rule combine_aggregated_weights:
    input: 
        expand(join(output_dir, 'weights/individual/{dataset}/gene_sets/{{agg_func}}/{feature}/{phenotype}.tsv.gz'),
               zip,
               #dataset=datasets, feature=features, phenotype=phenotypes)
        dataset=agg_datasets, feature=agg_features, phenotype=agg_phenotypes)
    output:
        join(output_dir, 'weights/combined/gene_sets/{agg_func}/combined_weights_{cor_method}_{collapse_func}.tsv.gz')
    script: 'src/combine/combine_weights.R'

rule build_individual_weights:
    input:
        features=join(output_dir, 'datasets/{dataset}/features/genes/{feature}.tsv.gz'),
        phenotype=join(output_dir, 'datasets/{dataset}/phenotypes/{phenotype}.tsv.gz')
    output:
        join(output_dir, 'weights/individual/{dataset}/genes/{feature}/{phenotype}.tsv.gz')
    script: 'src/datasource/pharmacogx/build_weights.R'

rule build_aggregated_weights:
    input:
        features=join(output_dir, 'datasets/{dataset}/features/gene_sets/{agg_func}/{feature}_gene_sets.tsv.gz'),
        phenotype=join(output_dir, 'datasets/{dataset}/phenotypes/{phenotype}.tsv.gz')
    output:
        join(output_dir, 'weights/individual/{dataset}/gene_sets/{agg_func}/{feature}/{phenotype}.tsv.gz')
    script: 'src/datasource/pharmacogx/build_weights.R'

rule merge_aggregated_features:
    input:
        expand(join(output_dir,
        'datasets/{{dataset}}/features/gene_sets/{gene_set}/{{agg_func}}/{{feature}}.tsv.gz'),
        gene_set=gene_set_names)
    output:
        join(output_dir, 'datasets/{dataset}/features/gene_sets/{agg_func}/{feature}_gene_sets.tsv.gz')

rule aggregate_features:
    input:
        join(output_dir, 'datasets/{dataset}/features/genes/{feature}.tsv.gz')
    output:
        join(output_dir, 'datasets/{dataset}/features/gene_sets/{gene_set}/{agg_func}/{feature}.tsv.gz')
    script: 'src/aggregate/aggregate_features.R'

rule parse_features:
    output:
        join(output_dir, 'datasets/{dataset}/features/genes/{feature}.tsv.gz')
    params:
        config=lambda wildcards, output: configs[wildcards.dataset]
    script: 'src/datasource/pharmacogx/parse_features.R'

rule parse_phenotypes:
    output:
        join(output_dir, 'datasets/{dataset}/phenotypes/{phenotype}.tsv.gz')
    params:
        config=lambda wildcards, output: configs[wildcards.dataset]
    script: 'src/datasource/pharmacogx/parse_phenotypes.R'

