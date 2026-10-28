import scanpy as sc
import hotspot
import pandas as pd
import matplotlib.pyplot as plt

import argparse
import time

parser = argparse.ArgumentParser(description='Process some files.')
parser.add_argument('--adata_file', type=str, required=True, help='Path to the adata file')
parser.add_argument('--output_folder', type=str, required=True, help='Path to the output directory')
parser.add_argument('--n_neighbors', type=int, default=8, help='n_neighbors')
parser.add_argument('--min_gene_threshold', type=int, default=50, help='min_gene_threshold')
parser.add_argument('--autoFDR', type=float, default=0.05, help='auto_FDR')
parser.add_argument(
    '--core_only',
    type=str,
    choices=['yes', 'no'],
    default='yes',
    help='only_keep core module yes / no'
)

args = parser.parse_args()

# 3. Convert yes/no to Boolean values
core_only_flag = True if args.core_only.lower() == 'yes' else False

# Load the Scanpy object
adata = sc.read_h5ad(args.adata_file)
sc.pp.calculate_qc_metrics(
    adata, percent_top=None, log1p=False, inplace=True
)
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=2)

#adata.obs['metabolite_counts'] = metabolite_counts[adata.obs.index.astype(int)]
adata.layers['raw'] = adata.X
adata.obs.head()
hs = hotspot.Hotspot(
    adata,
    layer_key='raw',
    model='danb',
    latent_obsm_key="spatial",
    umi_counts_obs_key="total_counts"
)
hs.create_knn_graph(weighted_graph=True, n_neighbors=args.n_neighbors)


hs_results = hs.compute_autocorrelations()
# Filter significant hotspot genes
hs_results.to_csv(f"{args.output_folder}/hs_result_danb_{args.n_neighbors}.csv")
#hs_genes = hs_results.index

local_correlations = hs.compute_local_correlations(adata.var_names, jobs=8) # jobs for parallelization

local_correlations.to_csv(f"{args.output_folder}/local_correlations_all_danb_{args.n_neighbors}.csv")
modules = hs.create_modules(
    min_gene_threshold=args.min_gene_threshold, core_only=core_only_flag, fdr_threshold=args.autoFDR # Set the significance threshold for merging two branches 
)
# Plot and save local correlations
plt.figure(figsize=(20, 20))  # For example, set the width to 12 inches and the height to 8 inches
hs.plot_local_correlations(vmin=-40,vmax=40)
plt.savefig(f"{args.output_folder}/hotmatrix_danb_n{args.n_neighbors}.pdf")


import joblib
joblib.dump(hs,f"{args.output_folder}/hs_danb.pkl")
