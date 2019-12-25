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


# wildcard_names = ['data_sources', 'feature_types', 'feature_levels']
#
# feature_wildcards = {wildcard: [] for wildcard in wildcards}
#
# for data_source_id in configs:
#     for feature_type in configs[data_source_id]['features']:
#         for feature_level in configs[data_source_id]['features'][feature_type]:
#             feature_wildcards['data_sources'].append(data_source_id)
#             feature_wildcards['feature_levels'].append(feature_level)
#             feature_wildcards['feature_types'].append(feature_type)

# weight combinations
cor_methods = config['combine']['correlation']['methods']
collapse_funcs = config['combine']['correlation']['collapse_funcs']

# create directory to store environment states (debug-mode)
if config['dev_mode']['enabled']:
    os.makedirs(config['dev_mode']['rda_dir'], mode = 755, exist_ok=True)

# rule all:
#     input:
#         expand(join(report_dir, '{feature_type}/{feature_level}/feature_weights_summary.html'),
#                feature_type=feature_types, feature_level=feature_levels)
#
# rule summarize_combined_weights:
#     input:
#         indiv_weights = expand(join(output_dir,
#             'weights/data_sources/{data_source}/{{feature_type}}/{{feature_level}}/{phenotype}.tsv.gz'), zip,
#             data_source=data_sources, phenotype=phenotypes),
#         combined_weights = expand(join(output_dir, 'weights/combined/{{feature_type}}/{{feature_level}}/combined_weights_{cor_method}_{collapse_func}.tsv.gz'),
#                                   cor_method=cor_methods, collapse_func=collapse_funcs)
#     output:
#         join(report_dir, '{feature_type}/{feature_level}/feature_weights_summary.html')
#     script:
#         'reports/feature_weights_summary.Rmd'
#
# rule combine_data_source_weights:
#     input:
#         expand(join(output_dir, 'weights/data_sources/{data_source}/{{feature_type}}/{{feature_level}}/{phenotype}.tsv.gz'),
#                zip,
#                data_source=data_sources, phenotype=phenotypes)
#     output:
#         join(output_dir, 'weights/combined/{feature_type}/{feature_level}/combined_weights_{cor_method}_{collapse_func}.tsv.gz')
#     script: 'src/combine/combine_weights.R'
#
# rule build_data_source_weights:
#     input:
#         features=join(output_dir, 'data_sources/{data_source}/features/{feature_type}/{feature_level}.tsv.gz'),
#         phenotype=join(output_dir, 'data_sources/{data_source}/phenotypes/{phenotype}.tsv.gz')
#     output:
#         join(output_dir, 'weights/data_sources/{data_source}/{feature_type}/{feature_level}/{phenotype}.tsv.gz')
#     script: 'src/weights/build_weights.R'

# rule all:
#     input:
#         expand(join(output_dir,
#               'weights/data_sources/{data_source}/{feature_type}/{feature_level}/{phenotype}/{cor_method}.feather'),
#                zip, feature_type=feature_types, feature_level=feature_levels)
#     output:
#         touch(join(output_dir, 'finished.txt'))

rule all:
    input:
        expand(join(output_dir, 'weights/combined/{feature_level}/combined_weights_{cor_method}_{collapse_func}.tsv.gz'), 
               feature_level=feature_levels, cor_method=cor_methods, collapse_func=collapse_funcs)

rule combine_data_source_weights:
    input:
        expand(join(output_dir, 'weights/data_sources/{data_source}/{feature_type}/{{feature_level}}/{phenotype}/{{cor_method}}_{{collapse_func}}.tsv.gz'),
               zip,
               data_source=data_sources, feature_type=feature_types, phenotype=phenotypes)
    output:
        join(output_dir, 'weights/combined/{feature_level}/combined_weights_{cor_method}_{collapse_func}.tsv.gz')
    script: 'src/combine/combine_weights.R'

rule compute_feature_phenotype_weights:
    input:
        cor_mat=join(output_dir, 'correlations/data_sources/{data_source}/{feature_type}/{feature_level}/{phenotype}/{cor_method}_cor_mat.feather')
    output:
        join(output_dir, 'weights/data_sources/{data_source}/{feature_type}/{feature_level}/{phenotype}/{cor_method}_{collapse_func}.tsv.gz')
    params:
        config=lambda wildcards, output: configs[wildcards.data_source]
    script: 'src/weights/compute_feature_weights.py'

rule compute_feature_phenotype_correlations:
    input:
        features=join(output_dir, 'data_sources/{data_source}/features/{feature_type}/{feature_level}.tsv.gz'),
        phenotype=join(output_dir, 'data_sources/{data_source}/phenotypes/{phenotype}.tsv.gz')
    output:
        join(output_dir, 'correlations/data_sources/{data_source}/{feature_type}/{feature_level}/{phenotype}/{cor_method}_cor_mat.feather')
    script: 'src/weights/compute_correlations.py'

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

