"""
" Compute feature-phenotype weights
"""
import numpy as np
import pandas as pd
import feather

# load saved correlation matrix
cor_mat = feather.read_dataframe(snakemake.input[0])

feat_type = snakemake.wildcards['feature_type']
feat_level = snakemake.wildcards['feature_level']

feat_cols = snakemake.params['config']['features'][feat_type][feat_level]['keys']

cor_mat = cor_mat.set_index(feat_cols)

# helper functions; will be moved to external file..
def maxabs(x):
    return max(abs(x))

# apply aggregation function along the correlation matrix rows
collapse_func = snakemake.wildcards['collapse_func']

if collapse_func == 'maxabs':
    res = cor_mat.apply(maxabs, axis=1)
else:
    res = cor_mat.apply(collapse_func, axis=1)

res = pd.DataFrame(res)
res.columns = [collapse_func]

res.to_csv(snakemake.output[0], sep='\t')

