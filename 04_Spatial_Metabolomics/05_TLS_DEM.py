# %%
import scanpy as sc
import anndata as ad
from collections import Counter
import pandas as pd
import gc
import matplotlib.pyplot as plt
import numpy as np
import os

# %%
lymphocyte_niche_sample_list = [
    "P36","P38","P32","P27","P39","P47","P46","P33","P49","P34","P37","P43_2","P45","P5","P8","P26",
    "P4","P43_1","P42","P40_2","P11","P30","P10"
]

# %% follicle-like TLS
ftls_sample_list = [
    "P38","P32","P27","P47","P46","P33","P49","P34","P37","P43_2"
]

#%% follicle-like TLS information
ftls_info = pd.read_csv("E:/Projects/BRCA/results/ST/Niche2607/TLS/TLS_cluster_classification.csv",
                        index_col=0)
ftls_info = ftls_info[ftls_info['maturation_cluster'] == 'fTLS']
ftls_clusters = ftls_info['cluster_full'].tolist()
ftls_clusters

# %% Compare follicle-like TLS with the remaining regions
for spID in ftls_sample_list[0:10]:
    print(spID)
    adata = sc.read_h5ad(f"E:/Projects/BRCA/Intermediate_Data/SM5u/matched_object/{spID}_match_st.h5ad")
    genes_delete = pd.read_csv('E:/Projects/BRCA/results/SM/filter/noise_metabolites_filtered.csv', index_col=0)
    genes_delete = genes_delete['gene'].astype(str).tolist()
    gene_mask = ~(adata.var_names.isin(genes_delete))
    adata = adata[:,gene_mask]
    
    sp_tls_cluster = pd.read_csv(f"E:/Projects/BRCA/results/ST/Niche2607/TLS/{spID}/TLS_dbscan_cluster_info.csv",
                             index_col=0)
    sp_tls_cluster.set_index('bin', inplace=True)
    sp_tls_cluster['dbscan_cluster_full'] = spID + "_" + sp_tls_cluster['dbscan_cluster']
    sp_tls_cluster['ftls_cluster'] = 'others'
    sp_tls_cluster.loc[sp_tls_cluster["dbscan_cluster_full"].isin(ftls_clusters), 'ftls_cluster'] = 'fTLS'
    sp_tls_cluster = sp_tls_cluster[sp_tls_cluster["ftls_cluster"] == "fTLS"]
    
    adata.obs['tls_cluster'] = 'others'
    adata.obs.loc[adata.obs.index.isin(sp_tls_cluster.index.tolist()),'tls_cluster'] = 'fTLS'
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.tl.rank_genes_groups(adata, "tls_cluster", method="wilcoxon")
    df  = sc.get.rank_genes_groups_df(adata, group="fTLS")
    df.to_csv(f"E:/Projects/BRCA/results/SM/TLS_2607/fTLS_dem_{spID}.csv")
    

# %% Compare follicle-like TLS with lymphocyte aggregate regions
for spID in ftls_sample_list[0:10]:
    print(spID)
    adata = sc.read_h5ad(f"E:/Projects/BRCA/Intermediate_Data/SM5u/matched_object/{spID}_match_st.h5ad")
    genes_delete = pd.read_csv('E:/Projects/BRCA/results/SM/filter/noise_metabolites_filtered.csv', index_col=0)
    genes_delete = genes_delete['gene'].astype(str).tolist()
    gene_mask = ~(adata.var_names.isin(genes_delete))
    adata = adata[:,gene_mask]
    
    sp_tls_cluster = pd.read_csv(f"E:/Projects/BRCA/results/ST/Niche2607/TLS/{spID}/TLS_dbscan_cluster_info.csv",
                             index_col=0)
    sp_tls_cluster.set_index('bin', inplace=True)
    sp_tls_cluster['dbscan_cluster_full'] = spID + "_" + sp_tls_cluster['dbscan_cluster']
    sp_tls_cluster['ftls_cluster'] = 'lymphocyte_region'
    sp_tls_cluster.loc[sp_tls_cluster["dbscan_cluster_full"].isin(ftls_clusters), 'ftls_cluster'] = 'fTLS'
    
    sp_ftls_cluster = sp_tls_cluster[sp_tls_cluster['ftls_cluster'] == 'fTLS']
    sp_lymphocyte_region = sp_tls_cluster[sp_tls_cluster['ftls_cluster'] == 'lymphocyte_region']
    
    if Counter(sp_tls_cluster['ftls_cluster'])['lymphocyte_region'] < 500:
        print("not enough lymphocyte region")
        continue
    
    adata.obs['tls_region'] = 'others'
    adata.obs.loc[adata.obs.index.isin(sp_ftls_cluster.index.tolist()),'tls_region'] = 'fTLS'
    adata.obs.loc[adata.obs.index.isin(sp_lymphocyte_region.index.tolist()),'tls_region'] = 'lymphocyte_region'
    print(Counter(adata.obs['tls_region']))
    
    
    adata = adata[adata.obs['tls_region'].isin(['fTLS','lymphocyte_region'])]
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.tl.rank_genes_groups(adata, "tls_region", method="wilcoxon")

    df1  = sc.get.rank_genes_groups_df(adata, group="fTLS")
    df1.to_csv(f"E:/Projects/BRCA/results/SM/fTLS_vs_lymphocyte_region_2607/fTLS_vs_lymphocyte_region_dem_{spID}.csv")
    
    df2 = sc.get.rank_genes_groups_df(adata, group="lymphocyte_region")
    df2.to_csv(f"E:/Projects/BRCA/results/SM/fTLS_vs_lymphocyte_region_2607/lymphocyte_region_vs_fTLS_dem_{spID}.csv")
# %%
