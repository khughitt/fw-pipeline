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

# add data source info to main config object
config['data_sources'] = configs

# Generate lists of output filepath components
data_sources = []
feature_types = []
feature_levels = []
phenotypes = []

for data_source_id in configs:
    for feature_type in configs[data_source_id]['features']:
        for feature_level in configs[data_source_id]['features'][feature_type]:
            for phenotype in configs[data_source_id]['phenotypes']:
                data_sources.append(data_source_id)
                feature_levels.append(feature_level)
                feature_types.append(feature_type)
                phenotypes.append(phenotype)

# weight combinations
cor_methods = config['combine']['correlation']['methods']
collapse_funcs = config['combine']['correlation']['collapse_funcs']

# create directory to store environment states (debug-mode)
if config['dev_mode']['enabled']:
    os.makedirs(config['dev_mode']['rda_dir'], mode = 755, exist_ok=True)

rule summarize_combined_weights:
    input: 
        indiv_weights = expand(join(output_dir,
            'weights/data_sources/{data_source}/{feature_type}/{feature_level}/{phenotype}.tsv.gz'), zip,
            data_source=data_sources, feature_type=feature_types, feature_level=feature_levels, phenotype=phenotypes),
        combined_weights = expand(join(output_dir, 'weights/combined/combined_weights_{cor_method}_{collapse_func}.tsv.gz'), cor_method=cor_methods, collapse_func=collapse_funcs),
    output:
        join(report_dir, 'combined_weights_genes_summary.html')
    script:
        'reports/combined_weights_summary.Rmd'

rule combine_data_source_weights:
    input: 
        expand(join(output_dir, 'weights/data_sources/{data_source}/{feature_type}/{feature_level}/{phenotype}.tsv.gz'),
               zip,
               data_source=data_sources, feature_type=feature_types, feature_level=feature_levels, phenotype=phenotypes)
    output:
        join(output_dir, 'weights/combined/combined_weights_{cor_method}_{collapse_func}.tsv.gz')
    params:
        feature_level='genes'
    script: 'src/combine/combine_weights.R'

rule build_data_source_weights:
    input:
        features=join(output_dir, 'data_sources/{data_source}/features/{feature_type}/{feature_level}.tsv.gz'),
        phenotype=join(output_dir, 'data_sources/{data_source}/phenotypes/{phenotype}.tsv.gz')
    output:
        join(output_dir, 'weights/data_sources/{data_source}/{feature_type}/{feature_level}/{phenotype}.tsv.gz')
    script: 'src/weights/build_weights.R'

rule load_features:
    output:
        join(output_dir, 'data_sources/{data_source}/features/{feature_type}/{feature_level}.tsv.gz')
    params:
        config=lambda wildcards, output: configs[wildcards.data_source]
    script: 'src/datasource/load_features.R'

rule load_phenotypes:
    output:
        join(output_dir, 'data_sources/{data_source}/phenotypes/{phenotype}.tsv.gz')
    params:
        config=lambda wildcards, output: configs[wildcards.data_source]
    script: 'src/datasource/load_phenotypes.R'

