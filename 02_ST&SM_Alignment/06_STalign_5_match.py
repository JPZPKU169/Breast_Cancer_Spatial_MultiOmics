import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
import json
import argparse

parser = argparse.ArgumentParser(description="Process input folder.")
parser.add_argument("--folder", type=str, help="Path to the input folder (positional)")
parser.add_argument("--spID", type=str, help="Sample ID (prefix of st cell id)")
args = parser.parse_args()

folder = args.folder
spID = args.spID

r_um = 2.5



# Load data
print("Loading data...")
df_spots = pd.read_csv(f"{folder}/3.3.target2source.csv")
df_cell_bin = pd.read_parquet(f"{folder}/barcode_mappings.parquet")
df_tissue = pd.read_parquet(f"{folder}/tissue_positions_2um.parquet")

# Preprocess: remove missing values and filter target cells
df_cell_bin_clean = df_cell_bin.dropna(subset=['cell_id']) # Remove 2 um spatial transcriptomics bins that are not assigned to a cell
df_cell_bin_clean['cell_id'] = spID + '_' + df_cell_bin_clean['cell_id']

# Preprocess: convert full-resolution image coordinates to high-resolution image coordinates
with open(f"{folder}/scalefactors_json.json", 'r') as f:
    data = json.load(f)

tissue_hires_scalef = data['tissue_hires_scalef']

df_tissue['pxl_row_in_hires'] = df_tissue['pxl_row_in_fullres'] * tissue_hires_scalef
df_tissue['pxl_col_in_hires'] = df_tissue['pxl_col_in_fullres'] * tissue_hires_scalef

# Merge to obtain bin coordinates
df_bin_coords = pd.merge(df_cell_bin_clean, df_tissue, left_on='square_002um', right_on='barcode', how='inner')

bin_coords = df_bin_coords[['pxl_col_in_hires', 'pxl_row_in_hires']].values
tree = cKDTree(bin_coords)

# Get spot coordinates
spot_coords = df_spots[['xJ_LDDMM', 'yJ_LDDMM']].values

# Query nearest neighbors
distances, indices = tree.query(spot_coords, k=1)

r = data['spot_diameter_fullres'] / 2 * tissue_hires_scalef * 2.5

matches = []
for i, (distance, idx) in enumerate(zip(distances, indices)):
    if distance < r:
        bin_id = df_bin_coords.iloc[idx]['square_002um']
        cell_id = df_bin_coords.iloc[idx]['cell_id']
        spot_id = df_spots.iloc[i]['cell']
        matches.append({'sm_id': spot_id, 'cell_id': cell_id})

# Create and save the results DataFrame
df_matches = pd.DataFrame(matches)
df_matches.to_csv(f"{folder}/sm_st_match.csv", index=False)
print(f"Matching complete. Results saved to {folder}/sm_st_match.csv")
