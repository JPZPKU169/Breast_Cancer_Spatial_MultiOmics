import argparse
import os
import scanpy as sc
from pool import microcluster, pool_matrix_cols
from scipy import sparse
import h5py
import pandas as pd
import numpy as np
from scipy import sparse

parser = argparse.ArgumentParser(description='micropool')
parser.add_argument('--expr_folder', type=str, required=True, help='Path to the adata folder')
args = parser.parse_args()

expr_folder = args.expr_folder

for fname in os.listdir(expr_folder):
    if fname.endswith('1.5.pseudo_expr.csv'):
        expr_path = os.path.join(expr_folder, fname)
        break
else:
    raise FileNotFoundError("No file ending with '1.5.pseudo_expr.csv' found in the folder.")

print(f"Loading expression data from: {expr_path}")

expr_df=pd.read_csv(expr_path, index_col=0)

adata= sc.AnnData(
    sparse.csr_matrix(expr_df.T.values),
    obs=pd.DataFrame(index=expr_df.columns), 
    var=pd.DataFrame(index=expr_df.index)
)

sc.settings.figdir=expr_folder
sc.pp.log1p(adata) 


sc.pp.highly_variable_genes(adata, n_top_genes=2000)

MIR = pd.read_csv("/n02dat02/PROJECT/liujianyu/ST/seurat/M_I_T.txt", sep="\t", header=0)
# Remove genes present in MIR['Gene.name'] from highly variable genes
print('MIR',MIR.columns.tolist())

mir_genes = set(MIR['Gene name'])
hv_genes = adata.var[adata.var['highly_variable']].index
filtered_hv_genes = [gene for gene in hv_genes if gene not in mir_genes]
# Update highly_variable column
adata.var['highly_variable'] = adata.var.index.isin(filtered_hv_genes)
print('Variable features after filtering:')
print(len(filtered_hv_genes))

sc.pl.highly_variable_genes(adata, save='_1.6.pdf')

if sparse.issparse(adata.X):
    adata.layers["scaled"] = adata.X.toarray().copy() 
else:
    adata.layers["scaled"] = adata.X.copy()

sc.pp.scale(adata, max_value=10, layer="scaled")

sc.pp.pca(adata,layer="scaled", svd_solver="arpack")

sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, save='_1.6.pdf')

sc.pl.pca(
    adata,
    dimensions=[(0, 1)],
    size=2,save='_1.6.pdf',
)


sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

sc.tl.umap(adata)

sc.tl.leiden(adata, resolution=1)

sc.pl.umap(adata, color=["leiden"],save='_1.6.pdf')

sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")

sc.pl.rank_genes_groups_dotplot(
    adata, groupby="leiden", standard_scale="var",min_logfoldchange=2 ,n_genes=7, save='1.6.pdf'
)

# Extract marker genes from rank_genes_groups results
markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.to_csv(f'{expr_folder}/1.6.all_markers.csv', index=True)

# For each cluster, select top 10 genes by log2 fold change
top7 = markers.groupby('group').head(7).reset_index(drop=True)
top7.to_csv(f'{expr_folder}/1.6.top_markers.csv', index=True)

adata.write(os.path.join(expr_folder, '1.6.pools_mapping.h5ad'))