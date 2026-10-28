import pandas as pd
import numpy as np
from scipy.spatial import KDTree
import argparse

def main():
    parser = argparse.ArgumentParser(description='Calculate cell type proportions with neighborhood filter')
    parser.add_argument('--celltype', required=True, help='Cell type CSV file')
    parser.add_argument('--coord', required=True, help='Coordinate CSV file')
    parser.add_argument('-r', '--radius', required=True, help='Radius value file (txt)')
    parser.add_argument('--output', required=True, help='Output CSV file')

    args = parser.parse_args()

    # Read the radius value
    with open(args.radius, 'r') as f:
        radius = float(f.read().strip())

    # Load data
    celltypes = pd.read_csv(args.celltype, index_col=0)
    coords = pd.read_csv(args.coord, index_col=0)
    
    # Get the list of shared bins
    common_bins = celltypes.index.intersection(coords.index)
    celltypes = celltypes.loc[common_bins]
    coords = coords.loc[common_bins]

    # Prepare data
    points = coords[['x', 'y']].values
    bin_ids = common_bins.values
    celltype_dict = celltypes['celltype'].to_dict()

    # Build the KDTree
    tree = KDTree(points)
    
    # Query all neighborhoods
    neighborhoods = tree.query_ball_tree(tree, radius)

    # Get all cell types
    all_types = celltypes['celltype'].unique().tolist()
    all_types.sort()

    # Calculate results
    results = []
    for i, neighbors in enumerate(neighborhoods):
        # Add a neighborhood size check
        #if len(neighbors) < 10:
        #    continue  # Skip this bin
        
        current_bin = bin_ids[i]
        
        # Get cell types in the neighborhood
        neighbor_types = [celltype_dict[bin_ids[j]] for j in neighbors]
        
        # Calculate proportions
        type_counts = pd.Series(neighbor_types).value_counts(normalize=True)
        
        # Construct the result row
        row = {'bin': current_bin}
        for ct in all_types:
            row[ct] = type_counts.get(ct, 0.0)
        
        results.append(row)

    # Create the results DataFrame, which may be empty
    if results:
        result_df = pd.DataFrame(results)
        result_df = result_df.reindex(columns=['bin'] + all_types)
        result_df.to_csv(args.output, index=False)
    else:
        pd.DataFrame(columns=['bin'] + all_types).to_csv(args.output, index=False)

if __name__ == '__main__':
    main()