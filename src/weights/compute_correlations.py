"""
" Compute feature-phenotype correlations
"""
from numba import njit
import numpy as np
import pandas as pd
import feather

"""
Efficient correlation implementation in python
Source: https://github.com/ikizhvatov/efficient-columnwise-correlation/blob/master/columnwise_corrcoef_perf.py
"""
def corrcoeff_einsum_optimized(O, P):
    (n, t) = O.shape      # n traces of t samples
    (n_bis, m) = P.shape  # n predictions for each of m candidates

    DO = O - (np.einsum("nt->t", O, optimize='optimal') / np.double(n)) # compute O - mean(O)
    DP = P - (np.einsum("nm->m", P, optimize='optimal') / np.double(n)) # compute P - mean(P)

    cov = np.einsum("nm,nt->mt", DP, DO, optimize='optimal')

    varP = np.einsum("nm,nm->m", DP, DP, optimize='optimal')
    varO = np.einsum("nt,nt->t", DO, DO, optimize='optimal')
    tmp = np.einsum("m,t->mt", varP, varO, optimize='optimal')

    return cov / np.sqrt(tmp)

"""
Fast Spearman correlation calculation

https://stackoverflow.com/questions/52371329/fast-spearman-correlation-between-two-pandas-dataframes
"""
@njit
def mean1(a):
  n = len(a)
  b = np.empty(n)
  for i in range(n):
    b[i] = a[i].mean()
  return b

@njit
def std1(a):
  n = len(a)
  b = np.empty(n)
  for i in range(n):
    b[i] = a[i].std()
  return b

@njit
def spearman_cor(a, b):
    """
    Spearman correlation

    Expects value rankings from unrolled 2d arrays with the same numbers of columns
    """
    n, k = a.shape
    m, k = b.shape

    mu_a = mean1(a)
    mu_b = mean1(b)
    sig_a = std1(a)
    sig_b = std1(b)

    out = np.empty((n, m))

    for i in range(n):
        for j in range(m):
            out[i, j] = (a[i] - mu_a[i]) @ (b[j] - mu_b[j]) / k / sig_a[i] / sig_b[j]

    return out

# main
#  mat = np.load(snakemake.input[0])

# store correlation matrix result for comparison
#  cor_mat = corrcoeff_einsum_optimized(mat.T, mat.T)
#  cor_mat[np.tril_indices_from(cor_mat)] = np.nan
#  np.save(snakemake.output['cor_mat'], cor_mat)

#corrcoeff_einsum_optimized(mat.T, mat.T)"

# get datasource settings
data_source_config = snakemake.config['data_sources'][snakemake.wildcards['data_source']]

# get correlation settings
correlation_params = snakemake.config['combine']['correlation']

# load processed feature / phenotype data
feat_dat = pd.read_csv(snakemake.input['features'], sep='\t')
pheno_dat = pd.read_csv(snakemake.input['phenotype'], sep='\t')

# determine feature sample indices
feat_type = snakemake.wildcards['feature_type']
feat_level = snakemake.wildcards['feature_level']

feat_keys = data_source_config['features'][feat_type][feat_level]['keys']

feat_sample_ids = [x for x in feat_dat.columns if x not in feat_keys] 

# match phenotype and feature sample ids
sample_ids = pheno_dat.columns[1:]
sample_ids = sorted(list(set(feat_sample_ids).intersection(sample_ids)))

pheno_key = pheno_dat.columns[0]

feat_dat = feat_dat.loc[:, feat_keys + sample_ids]
pheno_dat = pheno_dat.loc[:, [pheno_key] + sample_ids]

# generate numpy versions of the datasets without the id columns and transpose
feat_mat = feat_dat.loc[:, sample_ids].to_numpy().T
pheno_mat = pheno_dat.loc[:, sample_ids].to_numpy().T

# pearson correlation
if snakemake.wildcards['cor_method'] == 'pearson':
    cor_mat = corrcoeff_einsum_optimized(feat_mat, pheno_mat)
elif snakemake.wildcards['cor_method'] == 'spearman':
    cor_mat = spearman_cor(pd.DataFrame(feat_mat.T).rank(1).values, 
                           pd.DataFrame(pheno_mat.T).rank(1).values)
    cor_mat = cor_mat.T
else:
    raise Exception("Invalid correlation method specified!")

# convert back to a dataframe and add row/column identifiers
cor_mat = pd.DataFrame(cor_mat.T, columns=pheno_dat.drug)
cor_mat = pd.concat([feat_dat.loc[:, feat_keys], cor_mat], axis=1)

# store result
feather.write_dataframe(cor_mat, snakemake.output[0])
