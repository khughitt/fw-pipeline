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

# weight combinations
cor_methods = config['combine']['correlation']['methods']
collapse_funcs = config['combine']['correlation']['collapse_funcs']

# create directory to store environment states (debug-mode)
if config['dev_mode']['enabled']:
    os.makedirs(config['dev_mode']['rda_dir'], mode = 755, exist_ok=True)

rule run_fgsea_combined:
    input: 
        #expand(join(config['output_dir'], 'weights/combined_weights_{cor_method}_{collapse_func}s.tsv.gz'), 
        #       cor_method=cor_methods, collapse_func=collapse_funcs)
        join(config['output_dir'], 'weights/combined_weights_{cor_method}_{collapse_func}.tsv.gz')
    output:
        join(config['output_dir'], 'fgsea/combined_weights_{cor_method}_{collapse_func}_fgsea.tsv.gz')
    script: 'src/enrichment/run_fgsea.R'

rule summarize_fgsea_results:
    input: 
        indiv_weights = expand(join(config['output_dir'], 'weights/{dataset}/{version}/genes/{feature}/{phenotype}.tsv.gz'),
                               zip,
                               dataset=datasets, version=versions, feature=features, phenotype=phenotypes),
        indiv_fgsea = expand(join(config['output_dir'], 'fgsea/{dataset}/{version}/genes/{feature}/{phenotype}.tsv.gz'),
                             zip,
                             dataset=datasets, version=versions, feature=features, phenotype=phenotypes),
        combined_weights = expand(join(config['output_dir'], 'weights/combined_weights_{cor_method}_{collapse_func}.tsv.gz'),
                                  cor_method=cor_methods, collapse_func=collapse_funcs),
        combined_fgsea = expand(join(config['output_dir'], 'fgsea/combined_weights_{cor_method}_{collapse_func}_fgsea.tsv.gz'), 
                                cor_method=cor_methods, collapse_func=collapse_funcs)
    output:
        join(config['report_dir'], config['version'], 'fgsea_results.html')
    script:
        'reports/fgsea_results.Rmd'

rule run_fgsea_indiv:
    input: 
        join(config['output_dir'], 'weights/{dataset}/{version}/genes/{feature}/{phenotype}.tsv.gz') 
    output:
        join(config['output_dir'], 'fgsea/{dataset}/{version}/genes/{feature}/{phenotype}.tsv.gz')
    script: 'src/enrichment/run_fgsea.R'

rule combine_weights:
    input: 
        expand(join(config['output_dir'], 'weights/{dataset}/{version}/genes/{feature}/{phenotype}.tsv.gz'),
               zip,
               dataset=datasets, version=versions, feature=features, phenotype=phenotypes)
    output:
        join(config['output_dir'], 'weights/combined_weights_{cor_method}_{collapse_func}.tsv.gz')
    script: 'src/combine/combine_weights.R'

rule build_weights:
    input:
        features=join(config['output_dir'], 'datasets/{dataset}/{version}/genes/features/{feature}.tsv.gz'),
        phenotype=join(config['output_dir'], 'datasets/{dataset}/{version}/genes/phenotypes/{phenotype}.tsv.gz')
    output:
        join(config['output_dir'], 'weights/{dataset}/{version}/genes/{feature}/{phenotype}.tsv.gz')
    script: 'src/datasource/pharmacogx/build_weights.R'

rule parse_features:
    output:
        join(config['output_dir'], 'datasets/{dataset}/{version}/genes/features/{feature}.tsv.gz')
    params:
        config=lambda wildcards, output: configs[wildcards.dataset]
    script: 'src/datasource/pharmacogx/parse_features.R'

rule parse_phenotypes:
    output:
        join(config['output_dir'], 'datasets/{dataset}/{version}/genes/phenotypes/{phenotype}.tsv.gz')
    params:
        config=lambda wildcards, output: configs[wildcards.dataset]
    script: 'src/datasource/pharmacogx/parse_phenotypes.R'

