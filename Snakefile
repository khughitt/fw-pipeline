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

rule all:
    input: 
        expand(join(config['output_dir'], '{dataset}/{version}/genes/weights/{feature}/{phenotype}.tsv.gz'),
               dataset=datasets, version=versions, feature=features, phenotype=phenotypes)

rule build_weights:
    input:
        join(config['output_dir'], '{dataset}/{version}/genes/features/{feature}.tsv.gz'),
        join(config['output_dir'], '{dataset}/{version}/genes/phenotypes/{phenotype}.tsv.gz')
    output: join(config['output_dir'], '{dataset}/{version}/genes/weights/{feature}/{phenotype}.tsv.gz')
    shell: 'touch {output}'

rule parse_features:
    output: join(config['output_dir'], '{dataset}/{version}/genes/features/{feature}.tsv.gz')
    params:
        config=lambda wildcards, output: configs[wildcards.dataset]
    script: 'src/data_sources/pharmacogx_parse_features.R'

rule parse_phenotypes:
    output: join(config['output_dir'], '{dataset}/{version}/genes/phenotypes/{phenotype}.tsv.gz')
    params:
        config=lambda wildcards, output: configs[wildcards.dataset]
    script: 'src/data_sources/pharmacogx_parse_phenotypes.R'

