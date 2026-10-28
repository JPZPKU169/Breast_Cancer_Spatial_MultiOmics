import pandas as pd 
import numpy as np
import json
import torch
import pickle
from STalign import STalign
from matplotlib import pyplot as plt
import argparse

# load.STcoord
parser = argparse.ArgumentParser(description="Process input folder.")
parser.add_argument("--folder", type=str, help="Path to the input folder (positional)")

args = parser.parse_args()
folder = args.folder

# Load the spatial metabolomics spot coordinates
df=pd.read_csv(f'{folder}/coords.csv')
xJ = np.array(df['x'])
yJ = np.array(df['y']) 

with open(f'{folder}/2.1.out.pkl', 'rb') as f:
    out=pickle.load(f)
with open(f'{folder}/rasterized_IJ.pkl', 'rb') as f:
    XJ, YJ, J,XI,YI,I,J,muA,muB= pickle.load(f) ## Landmarks were initially selected on the CPU; move variables such as XI and YI to the GPU before subsequent processing
extentI = STalign.extent_from_x((YI,XI))

A_post = out['A']
v = out['v']
xv = out['xv']
WM = out['WM'] 
WB=out['WB']
WA=out['WA']
#J = torch.tensor(J, dtype=torch.float32, device=v.device)

extentJ = STalign.extent_from_x((YJ,XJ))

hirecoord=torch.tensor(np.stack([yJ,xJ], axis=1), dtype=torch.float32, device=v.device) # Use y before x


A_post = A_post.float()   # Convert A to float32
tpointsI= STalign.transform_points_target_to_source(xv,v,A_post,hirecoord)

if tpointsI.is_cuda:
    tpointsI = tpointsI.cpu()

padding_width = 300

xJ_LDDMM = tpointsI[:,1] - padding_width
yJ_LDDMM = tpointsI[:,0] - padding_width

# plot results
fig,ax = plt.subplots(figsize=(20, 20))
ax.imshow((I).transpose(1,2,0),extent=extentI)
#ax.imshow(phiI.permute(1,2,0)/torch.max(phiI))
ax.scatter(xJ_LDDMM,yJ_LDDMM,s=0.01,alpha=0.3,label = 'source aligned')
ax.set_aspect('equal')
lgnd = plt.legend(loc="upper right", scatterpoints=1, fontsize=10)
for handle in lgnd.legend_handles:
    handle.set_sizes([10.0])

ax.set_title('After alignment source and spot of source ', fontsize=15)
plt.savefig(f'{folder}/3.2.target2source.png',dpi=600)



# Combine locs with xI_LDDMM and yI_LDDMM
df['xJ_LDDMM'] = xJ_LDDMM.numpy()
df['yJ_LDDMM'] = yJ_LDDMM.numpy()
# Save the combined DataFrame to a CSV file
df.to_csv(f'{folder}/3.3.target2source.csv', index=False)

