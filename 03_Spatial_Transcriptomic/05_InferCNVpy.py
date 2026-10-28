import scanpy as sc
import infercnvpy as cnv
import matplotlib.pyplot as plt
import warnings
import pandas as pd
import argparse
import os
import numpy as np

parser = argparse.ArgumentParser(description="InferCNVpy")
parser.add_argument("--data_dir", type=str, default="/home/output", help="Output directory")
args = parser.parse_args()

data_dir=args.data_dir
output_dir=data_dir

sc.settings.set_figure_params(figsize=(5, 5))
sc.settings.figdir=output_dir

adata_epi_t = sc.read_h5ad(f"{data_dir}/metacell_infercnv_log1p_input.h5ad")

cnv.tl.infercnv(
    adata_epi_t,
    reference_key="celltype",
    reference_cat=[
        "T"
    ],
    window_size=250,
    chunksize=3000,
    n_jobs=1
)

cnv.pl.chromosome_heatmap(adata_epi_t, groupby="celltype", save=f"2.1.chromosome_heatmap.png",show=False)

###Clustering by CNV profiles and identifying tumor cells
cnv.tl.pca(adata_epi_t)
cnv.pp.neighbors(adata_epi_t)
cnv.tl.leiden(adata_epi_t)

cnv.pl.chromosome_heatmap(adata_epi_t, groupby="cnv_leiden", dendrogram=True)

cnv.tl.umap(adata_epi_t)
cnv.tl.cnv_score(adata_epi_t)


adata_epi_t.write_h5ad(f'{output_dir}/2.4.adata_infercnv.h5ad')

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 11))
ax4.axis("off")
cnv.pl.umap(
    adata_epi_t,
    color="cnv_leiden",
    legend_loc="on data",
    legend_fontoutline=2,
    ax=ax1,
    show=False,
)
cnv.pl.umap(adata_epi_t, color="cnv_score", ax=ax2, show=False)
cnv.pl.umap(adata_epi_t, color="celltype", ax=ax3)
fig.savefig(f'{output_dir}/2.2.cnv_umap_leiden.png', dpi=900)

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 11), gridspec_kw={"wspace": 0.5})
ax4.axis("off")
sc.pl.umap(adata_epi_t, color="cnv_leiden", ax=ax1, show=False)
sc.pl.umap(adata_epi_t, color="cnv_score", ax=ax2, show=False)
sc.pl.umap(adata_epi_t, color="celltype", ax=ax3)
fig.savefig(f'{output_dir}/2.3.sample_umap_leiden.png', dpi=900)



