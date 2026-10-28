import scanpy as sc
import pandas as pd
import numpy as np
from scipy.sparse import issparse, csr_matrix, lil_matrix
import anndata as ad
import argparse
from tqdm import tqdm

parser = argparse.ArgumentParser(description="Process input folder.")
parser.add_argument("--folder", type=str, help="Path to the input folder (positional)")
parser.add_argument("--spID", type=str, help="sample ID as prefix of ST cell ID")
parser.add_argument("--sm_adata", type=str, help="path to the h5ad file of spatial metabolomics")
parser.add_argument("--output_folder", type=str, help="path to save st matched sm object.")


args = parser.parse_args()
folder = args.folder
spID = args.spID

# Load the mapping between spatial metabolomics spots and spatial transcriptomics cells
df_matches = pd.read_csv(f"{folder}/sm_st_match.csv", dtype=str)

df_matches['sm_id'] = df_matches['sm_id'].str.replace(r'\.0$', '', regex=True) # Read as strings and remove the trailing ".0" from numeric IDs

cell_seg_position = pd.read_csv(f"{folder}/tissue_positions_seg.csv")

cell_seg_position['barcode'] = spID + "_" + cell_seg_position['barcode']

adata = sc.read_h5ad(args.sm_adata)

if issparse(adata.X):
    print("✓ The current matrix is already sparse")
    print(f"  Sparse matrix format: {adata.X.format}")
    print(f"  Number of nonzero elements: {adata.X.nnz}")
    print(f"  Matrix shape: {adata.X.shape}")
else:
    print("✗ The current matrix is dense; converting to a CSR sparse matrix...")
        
    # Convert to a CSR sparse matrix
    adata.X = csr_matrix(adata.X)

def efficient_merge_spots_to_cells(adata, df_matches, cell_seg_position, chunk_size=10000):
    """
    Efficiently aggregate large-scale spatial metabolomics data at the cell level
    
    Parameters:
    adata: AnnData object containing spatial metabolomics data
    df_matches: DataFrame containing sm_id-to-cell_id matches
    cell_seg_position: DataFrame containing spatial coordinates of cells
    chunk_size: Number of items per processing chunk
    
    Returns:
    An aggregated AnnData object
    """
    # 1. Filter adata to retain only matched spots
    matched_sm_ids = set(df_matches['sm_id'].tolist())
    mask = [name in matched_sm_ids for name in adata.obs_names]
    adata_filtered = adata[mask].copy()
    
    # 2. Add cell_id information to adata.obs
    sm_to_cell = dict(zip(df_matches['sm_id'], df_matches['cell_id']))
    adata_filtered.obs['cell_id'] = [sm_to_cell.get(name, 'unknown') for name in adata_filtered.obs_names]
    
    # 3. Get all unique cell_ids
    unique_cell_ids = sorted(adata_filtered.obs['cell_id'].unique())
    n_cells = len(unique_cell_ids)
    n_features = adata_filtered.shape[1]
    
    # 4. Create a mapping from cell_id to index
    cell_to_idx = {cell_id: idx for idx, cell_id in enumerate(unique_cell_ids)}
    
    # 5. Preallocate the result as a sparse matrix
    if hasattr(adata_filtered.X, 'todense'):
        result_matrix = lil_matrix((n_cells, n_features), dtype=adata_filtered.X.dtype)
    else:
        result_matrix = np.zeros((n_cells, n_features), dtype=adata_filtered.X.dtype)
    
    # 6. Group by cell_id and process in batches
    cell_groups = adata_filtered.obs.groupby('cell_id').groups
    
    for i, cell_id in enumerate(tqdm(unique_cell_ids, desc="Aggregating cell expression")):
        spot_indices = cell_groups[cell_id]
        
        cell_data = adata_filtered[spot_indices]
        
        if hasattr(cell_data.X, 'sum'):
            summed = np.array(cell_data.X.sum(axis=0)).flatten()
        else:
            summed = cell_data.X.sum(axis=0)
        
        cell_idx = cell_to_idx[cell_id]
        result_matrix[cell_idx] = summed
    
    # 7. Convert to the appropriate matrix format
    if hasattr(adata_filtered.X, 'todense'):
        result_matrix = result_matrix.tocsr()
    
    # 8. Create a new AnnData object
    new_adata = ad.AnnData(
        X=result_matrix,
        obs=pd.DataFrame(index=unique_cell_ids),
        var=adata_filtered.var.copy()
    )
    
    # 9. Add spatial coordinates
    cell_to_coords = {}
    for _, row in cell_seg_position.iterrows():
        cell_id = row['barcode']
        x_coord = row['pxl_col_in_fullres']
        y_coord = row['pxl_row_in_fullres']
        cell_to_coords[cell_id] = (x_coord, y_coord)
    
    new_adata.obs['x'] = new_adata.obs_names.map(lambda x: cell_to_coords.get(x, (np.nan, np.nan))[0])
    new_adata.obs['y'] = new_adata.obs_names.map(lambda x: cell_to_coords.get(x, (np.nan, np.nan))[1])
    
    # 10. Copy any other required attributes
    new_adata.uns = adata_filtered.uns.copy()
    
    return new_adata

adata_merged = efficient_merge_spots_to_cells(adata, df_matches, cell_seg_position, chunk_size=10000)
sc.pp.calculate_qc_metrics(
    adata_merged, percent_top=None, log1p=False, inplace=True
)

adata_merged.obsm['spatial'] = adata_merged.obs[['x','y']].values

smID = spID.replace("ST","SM")
sc.write(f"{args.output_folder}/{smID}_match_st.h5ad", adata_merged)