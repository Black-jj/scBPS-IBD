from pathlib import Path
import os
PACKAGE_ROOT = Path(__file__).resolve().parents[1]
INT = str(PACKAGE_ROOT / 'work' / '20_geneformer')
OUT = str(PACKAGE_ROOT / 'results' / '20_Geneformer')
GFSRC = str(PACKAGE_ROOT / 'data' / 'reference' / 'geneformer')
MODEL = str(PACKAGE_ROOT / 'data' / 'models' / 'Geneformer-V2-104M')
GENES = ['CCL20', 'CLDN4', 'CXCL8']
os.environ.setdefault('TRANSFORMERS_OFFLINE', '1')
for directory in (Path(INT) / 'h5ad', Path(INT) / 'gf_tokenized', Path(INT) / 'gf_perturb', Path(INT) / 'gf_stats', Path(OUT)):
    directory.mkdir(parents=True, exist_ok=True)
import os, pickle, pandas as pd
genes = [x.strip() for x in open(os.path.join(INT, 'gf_input', 'genes.txt'))]
sym2ens = pickle.load(open(os.path.join(GFSRC, 'gene_name_id_dict_gc104M.pkl'), 'rb'))
tokdict = pickle.load(open(os.path.join(GFSRC, 'token_dictionary_gc104M.pkl'), 'rb'))
rows = []
for symbol in genes:
    ens = sym2ens.get(symbol.upper())
    if ens in tokdict:
        rows.append((symbol, symbol.upper(), ens, tokdict[ens]))
mapping = pd.DataFrame(rows, columns=['input_gene_symbol', 'geneformer_gene_symbol', 'geneformer_ensembl_id', 'geneformer_token'])
required = set(GENES)
observed = set(mapping['geneformer_gene_symbol'])
if not required.issubset(observed):
    raise RuntimeError('Missing KO genes in Geneformer mapping: ' + ','.join(sorted(required - observed)))
mapping.to_csv(os.path.join(INT, 'gf_input', 'human_geneformer_mapping.tsv'), sep='\t', index=False)
import os, pickle, numpy as np, pandas as pd, scipy.io, scipy.sparse as sp, anndata as ad
M = scipy.io.mmread(os.path.join(INT, 'gf_input', 'matrix.mtx')).tocsr()
genes = [l.strip() for l in open(os.path.join(INT, 'gf_input', 'genes.txt'))]
barc = [l.strip() for l in open(os.path.join(INT, 'gf_input', 'barcodes.txt'))]
meta = pd.read_csv(os.path.join(INT, 'gf_input', 'meta.csv'))
X = M.T.tocsr()
sym2ens = pickle.load(open(os.path.join(GFSRC, 'gene_name_id_dict_gc104M.pkl'), 'rb'))
ens = [sym2ens.get(g.upper()) for g in genes]
keep = [i for i, e in enumerate(ens) if e is not None]
X = X[:, keep]
v_ens = [ens[i] for i in keep]
v_sym = [genes[i] for i in keep]
seen = set()
uniq = [i for i, e in enumerate(v_ens) if not (e in seen or seen.add(e))]
X = X[:, uniq]
v_ens = [v_ens[i] for i in uniq]
v_sym = [v_sym[i] for i in uniq]
required = set(GENES)
observed = set((x.upper() for x in v_sym))
if not required.issubset(observed):
    raise RuntimeError('Missing KO genes after human Geneformer mapping: ' + ','.join(sorted(required - observed)))
adata = ad.AnnData(X=X)
adata.var['ensembl_id'] = v_ens
adata.var['gene_symbol'] = v_sym
adata.var_names = v_sym
adata.obs['n_counts'] = np.asarray(X.sum(1)).ravel()
adata.obs['barcode'] = barc
for c in ['group3', 'group', 'cellType_1']:
    adata.obs[c] = meta[c].values
adata.write(os.path.join(INT, 'h5ad', 'enterocytes.h5ad'))
from geneformer import TranscriptomeTokenizer
tk = TranscriptomeTokenizer(
    {'barcode': 'barcode', 'group3': 'group3', 'group': 'group', 'cellType_1': 'cellType_1'},
    nproc=4,
    model_version='V2',
    gene_median_file=os.path.join(GFSRC, 'gene_median_dictionary_gc104M.pkl'),
    token_dictionary_file=os.path.join(GFSRC, 'token_dictionary_gc104M.pkl'),
    gene_mapping_file=os.path.join(GFSRC, 'ensembl_mapping_dict_gc104M.pkl')
)
tk.tokenize_data(os.path.join(INT, 'h5ad'), os.path.join(INT, 'gf_tokenized'), 'enterocytes', file_format='h5ad')
import os, pickle, glob
genes = GENES
sym2ens = pickle.load(open(os.path.join(GFSRC, 'gene_name_id_dict_gc104M.pkl'), 'rb'))
tokdict = pickle.load(open(os.path.join(GFSRC, 'token_dictionary_gc104M.pkl'), 'rb'))
from geneformer import InSilicoPerturber
ds = os.path.join(INT, 'gf_tokenized', 'enterocytes.dataset')
for sym in genes:
    lookup = sym.upper()
    ens = sym2ens[lookup]
    tok = tokdict[ens]
    outdir = os.path.join(INT, 'gf_perturb', sym)
    os.makedirs(outdir, exist_ok=True)
    done = glob.glob(os.path.join(outdir, f'*_gene_embs_dict_[[]{tok}[]]_raw.pickle'))
    if done:
        continue
    isp = InSilicoPerturber(
        perturb_type='delete', genes_to_perturb=[ens], model_type='Pretrained',
        emb_mode='cls_and_gene', max_ncells=None, emb_layer=-1,
        forward_batch_size=24, nproc=1, model_version='V2',
        token_dictionary_file=os.path.join(GFSRC, 'token_dictionary_gc104M.pkl')
    )
    isp.perturb_data(MODEL, ds, outdir, sym)
import os, pickle, glob, numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats as st
from statsmodels.stats.multitest import multipletests
genes = GENES
sym2ens = pickle.load(open(os.path.join(GFSRC, 'gene_name_id_dict_gc104M.pkl'), 'rb'))
tokdict = pickle.load(open(os.path.join(GFSRC, 'token_dictionary_gc104M.pkl'), 'rb'))
TOK = os.path.join(GFSRC, 'token_dictionary_gc104M.pkl')
NID = os.path.join(GFSRC, 'gene_name_id_dict_gc104M.pkl')
ligands = set((l.strip().upper() for l in open(os.path.join(INT, 'gf_input', 'ligands.txt'))))
tok2ens = {v: k for k, v in tokdict.items()}
ens2sym = {v: k for k, v in sym2ens.items()}
MINN = 50
statdir = os.path.join(INT, 'gf_stats')
pdf = PdfPages(os.path.join(OUT, '1.Geneformer_Volcano_All_KO.pdf'))
top_changed = []
ligand_rows = []
for sym in genes:
    lookup = sym.upper()
    ens = sym2ens[lookup]
    tok = tokdict[ens]
    pk = os.path.join(INT, 'gf_perturb', sym, f'in_silico_delete_{sym}_gene_embs_dict_[{tok}]_raw.pickle')
    d = pickle.load(open(pk, 'rb'))
    rows = []
    for (a, t), vals in d.items():
        if t == tok or len(vals) == 0:
            continue
        arr = np.asarray(vals, dtype=float)
        n = arr.size
        rows.append((t, float(arr.mean()), float(arr.std(ddof=1)) if n > 1 else 0.0, n))
    df = pd.DataFrame(rows, columns=['tok', 'Cosine_sim_mean', 'Cosine_sim_stdev', 'N_Detections'])
    df['Affected_Ensembl_ID'] = df['tok'].map(tok2ens)
    df['Geneformer_human_gene_name'] = df['Affected_Ensembl_ID'].map(ens2sym)
    df['Affected_gene_name'] = df['Geneformer_human_gene_name']
    df = df[df['Affected_gene_name'].notna()].copy()
    df = df[df['N_Detections'] >= MINN].copy()
    df = df.reset_index(drop=True)
    mu0 = np.median(df['Cosine_sim_mean'])
    se = df['Cosine_sim_stdev'] / np.sqrt(df['N_Detections'].clip(lower=1))
    se = se.replace(0, np.nan)
    tval = (df['Cosine_sim_mean'] - mu0) / se
    dfree = (df['N_Detections'] - 1).clip(lower=1)
    pval = 2 * st.t.sf(np.abs(tval), dfree)
    pval = np.where(np.isnan(pval), 1.0, pval)
    fdr = multipletests(pval, method='fdr_bh')[1]
    df['Cosine_Shift'] = mu0 - df['Cosine_sim_mean']
    df['pval'] = pval
    df['FDR'] = fdr
    df['neglog10FDR'] = -np.log10(np.clip(df['FDR'], 1e-300, None))
    df['significant'] = df['FDR'] < 0.05
    df = df.sort_values('Cosine_Shift', ascending=False).reset_index(drop=True)
    df.to_csv(os.path.join(statdir, f'{sym}_affected_stats.csv'), index=False)
    sig = df[df['significant']]
    if len(sig):
        top = sig.iloc[0]
        top_changed.append(dict(KO_gene=sym, Top_changed_gene=top['Affected_gene_name'], Cosine_Shift=round(float(top['Cosine_Shift']), 4), FDR=float(top['FDR']), N_Detections=int(top['N_Detections']), n_significant=int(len(sig))))
    else:
        top_changed.append(dict(KO_gene=sym, Top_changed_gene='NA', Cosine_Shift=np.nan, FDR=np.nan, N_Detections=0, n_significant=0))
    lig = sig[sig['Affected_gene_name'].str.upper().isin(ligands)].copy()
    lig['ligand_activity'] = lig['Cosine_Shift'] * lig['neglog10FDR']
    for _, r in lig.iterrows():
        ligand_rows.append(dict(KO_gene=sym, ligand=r['Affected_gene_name'], Cosine_Shift=round(float(r['Cosine_Shift']), 4), FDR=float(r['FDR']), neglog10FDR=round(float(r['neglog10FDR']), 3), ligand_activity=round(float(r['ligand_activity']), 4)))
    fig, ax = plt.subplots(figsize=(6, 5.2))
    ax.scatter(df['Cosine_Shift'], df['neglog10FDR'], s=10, c='#BBBBBB', edgecolors='none', label='ns')
    red = df[df['significant']]
    ax.scatter(red['Cosine_Shift'], red['neglog10FDR'], s=14, c='#D62728', edgecolors='none', label='FDR<0.05')
    ax.axhline(-np.log10(0.05), ls='--', lw=0.7, c='grey')
    ax.axvline(0, ls='--', lw=0.7, c='grey')
    xmin = float(df['Cosine_Shift'].min())
    xmax = float(df['Cosine_Shift'].max())
    xspan = max(xmax - xmin, 1e-06)
    ax.set_xlim(xmin - 0.04 * xspan, xmax + 0.14 * xspan)
    right_cutoff = float(df['Cosine_Shift'].quantile(0.9))
    for _, r in red.sort_values('Cosine_Shift', ascending=False).head(10).iterrows():
        align_right = r['Cosine_Shift'] >= right_cutoff
        ax.annotate(r['Affected_gene_name'], (r['Cosine_Shift'], r['neglog10FDR']), fontsize=7, xytext=(-3 if align_right else 3, 3), textcoords='offset points', ha='right' if align_right else 'left')
    ax.set_xlabel('Cosine Shift  (>0 = stronger shift than background)')
    ax.set_ylabel('$-\\log_{10}$ FDR')
    ax.set_title(f'Enterocytes  KO:{sym}  (affected genes, n={len(df)})')
    ax.legend(fontsize=7, loc='upper left', frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(statdir, f'volcano_{sym}.png'), dpi=200)
    fig.savefig(os.path.join(OUT, f'2.Geneformer_Volcano_KO_{sym}.pdf'))
    pdf.savefig(fig)
    plt.close(fig)
pdf.close()
pd.DataFrame(top_changed).to_csv(os.path.join(statdir, 'top_changed_gene_per_KO.csv'), index=False)
lig_df = pd.DataFrame(ligand_rows)
if len(lig_df):
    lig_df = lig_df.sort_values(['ligand_activity'], ascending=False).reset_index(drop=True)
lig_df.to_csv(os.path.join(statdir, 'top_ligand_activity.csv'), index=False)
