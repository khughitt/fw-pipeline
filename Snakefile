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

rule summarize_fgsea_results:
    input: 
        expand(join(config['output_dir'], '{dataset}/{version}/genes/fgsea/{feature}/{phenotype}.tsv.gz'),
               zip,
               dataset=datasets, version=versions, feature=features, phenotype=phenotypes),
        join(config['output_dir'], 'fgsea/combined_weights_fgsea.tsv.gz')
    output:
        join(config['report_dir'], config['version'], 'fgsea_results.html')
    script:
        'reports/fgsea_results.Rmd'

rule run_fgsea_combined:
    input: 
        join(config['output_dir'], 'weights/combined_weights.tsv.gz')
    output:
        join(config['output_dir'], 'fgsea/combined_weights_fgsea.tsv.gz')
    script: 'src/enrichment/run_fgsea.R'

rule run_fgsea_indiv:
    input: 
        join(config['output_dir'], '{dataset}/{version}/genes/weights/{feature}/{phenotype}.tsv.gz') 
    output:
        join(config['output_dir'], '{dataset}/{version}/genes/fgsea/{feature}/{phenotype}.tsv.gz')
    script: 'src/enrichment/run_fgsea.R'

rule combine_weights:
    input: 
        expand(join(config['output_dir'], '{dataset}/{version}/genes/weights/{feature}/{phenotype}.tsv.gz'),
               zip,
               dataset=datasets, version=versions, feature=features, phenotype=phenotypes)
    output:
        join(config['output_dir'], 'weights/combined_weights.tsv.gz')
    script: 'src/combine/combine_weights.R'

rule build_weights:
    input:
        features=join(config['output_dir'], '{dataset}/{version}/genes/features/{feature}.tsv.gz'),
        phenotype=join(config['output_dir'], '{dataset}/{version}/genes/phenotypes/{phenotype}.tsv.gz')
    output:
        join(config['output_dir'], '{dataset}/{version}/genes/weights/{feature}/{phenotype}.tsv.gz')
    script: 'src/datasource/pharmacogx/compute_cor.R'

rule parse_features:
    output:
        join(config['output_dir'], '{dataset}/{version}/genes/features/{feature}.tsv.gz')
    params:
        config=lambda wildcards, output: configs[wildcards.dataset]
    script: 'src/datasource/pharmacogx/parse_features.R'

rule parse_phenotypes:
    output:
        join(config['output_dir'], '{dataset}/{version}/genes/phenotypes/{phenotype}.tsv.gz')
    params:
        config=lambda wildcards, output: configs[wildcards.dataset]
    script: 'src/datasource/pharmacogx/parse_phenotypes.R'

