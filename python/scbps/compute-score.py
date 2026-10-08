import scdrs
from anndata import AnnData
from scipy import stats
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import warnings
import scanpy as sc
import gc
import joblib
import argparse
import logging
parser = argparse.ArgumentParser(description='compute-score')
parser.add_argument('input_file', type=str, help='Path to the input file (tab-separated).')
parser.add_argument('out_folder', type=str, help='Path to the output folder.')
parser.add_argument('adata_file', type=str, help='Path to the adata preprocessed file.')
args = parser.parse_args()
logfile = args.out_folder + 'compute_score.log'
logging.basicConfig(filename=logfile, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
def log_message(message):
    logging.info(message)
ff1 = ['%s', '%.5f', '%.8f', '%.4e', '%.4e', '%.5f', '%.6f']
ff2 = ['%.5f'] * 1000
ff = ff1 + ff2
def df2csv(df, fname, myformats=[], sep='\t'):
    if len(df.columns) <= 0:
        return
    Nd = len(df.columns)
    Nd_1 = Nd - 1
    formats = myformats[:]
    Nf = len(formats)
    if Nf < Nd:
        for ii in range(Nf, Nd):
            coltype = df[df.columns[ii]].dtype
            ff = '%s'
            if coltype == np.int64:
                ff = '%d'
            elif coltype == np.float64:
                ff = '%f'
            formats.append(ff)
    fh = open(fname, 'w')
    fh.write('\t'.join(df.columns) + '\n')
    for row in df.itertuples(index=False):
        ss = ''
        for ii in range(Nd):
            ss += formats[ii] % row[ii]
            if ii < Nd_1:
                ss += sep
        fh.write(ss + '\n')
    fh.close()
OUT_FOLDER = args.out_folder
if not os.path.exists(OUT_FOLDER):
    os.makedirs(OUT_FOLDER)
df_gs = scdrs.util.load_gs(args.input_file)
ls = pd.read_csv(args.input_file, sep='\t', index_col=0)
adata_file_path = args.adata_file
adata = sc.read_h5ad(adata_file_path)
scdrs.preprocess(adata)
for i in range(len(ls)):
    trait = ls.index[i]
    gene_list = df_gs[trait][0]
    gene_weight = df_gs[trait][1]
    df_res = scdrs.score_cell(adata, gene_list, gene_weight=gene_weight, n_ctrl=1000, return_ctrl_norm_score=False, copy=True)
    df_res.insert(0, 'cell_id', df_res.index)
    df2csv(df_res.iloc[:, 0:7], os.path.join(OUT_FOLDER, '%s.score' % trait), myformats=ff1)
    del df_res
    gc.collect()
log_message('Compute score finished')
