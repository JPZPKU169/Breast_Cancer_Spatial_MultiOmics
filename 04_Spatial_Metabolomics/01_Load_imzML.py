import pandas as pd
import numpy as np
from pyimzml.ImzMLParser import ImzMLParser
import anndata
import argparse
# Read the imzML file and specify the output location
parser = argparse.ArgumentParser(description='Process some files.')
parser.add_argument('--imzml_file', type=str, required=True, help='Path to the imzml file')
parser.add_argument('--output_file', type=str, required=True, help='Path to the output directory')

args = parser.parse_args()



#---------------------------------------------------

p = ImzMLParser(args.imzml_file)

# Extract pixel data
image_data = []
for idx, (x, y, z) in enumerate(p.coordinates):
    image_data.append({'x': x, 'y': y, 'z': z})

image_data = pd.DataFrame(image_data) 
image_data['Breaks'] = 1

count = 1
for i in range(1, len(image_data)):
    if abs(image_data['y'].iloc[i] - image_data['y'].iloc[i-1]) >= 2:
        count += 1
    image_data['Breaks'].iloc[i] = count
image_data['Breaks']=1

# Extract m/z and intensity data
mz_values = []

intensity = []

# Extract all spectra and construct the set of m/z values
for i in range(len(p.coordinates)):
    mzs, intensities = p.getspectrum(i)
    mz_values.append(mzs)
    intensity.append(intensities)

# Merge and sort all m/z values
all_mz_values = sorted(set(np.concatenate(mz_values)))

inten_matrix = pd.DataFrame(np.nan, index=range(len(mz_values)), columns=all_mz_values)

for i in range(len(mz_values)):
    inten_matrix.loc[i, mz_values[i]] = intensity[i]

inten_matrix = inten_matrix.fillna(0)

print(inten_matrix.head())
adata=[]
adata = anndata.AnnData(X=inten_matrix)
adata.var['mz'] = inten_matrix.columns
adata.obs['Breaks'] = pd.Categorical(image_data['Breaks'])

adata.obsm['spatial'] = image_data[['x', 'y']].values
adata.write(args.output_file)