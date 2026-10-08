import os
import pandas as pd
import argparse
parser = argparse.ArgumentParser(description='generate_norm_score')
parser.add_argument('infile', type=str, help='Path to the input file (tab-separated).')
parser.add_argument('outfile', type=str, help='Path to the output folder.')
args = parser.parse_args()
df_gs = pd.read_csv(args.infile, header=None)
dict_df_stats = {trait: pd.read_csv(trait, sep='\t', index_col=0) for trait in df_gs.iloc[:, 0]}
trait_list = list(dict_df_stats.keys())
df_norm = pd.concat([dict_df_stats[trait]['norm_score'] for trait in trait_list], axis=1)
df_norm.columns = trait_list
df_norm.columns = [os.path.basename(col).replace('.score', '') for col in df_norm.columns]
df_norm.to_csv(args.outfile, sep='\t')
