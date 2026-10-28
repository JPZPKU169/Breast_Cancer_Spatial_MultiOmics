import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import torch
from STalign import STalign
import pickle
import matplotlib.transforms as mtransforms
import argparse

parser = argparse.ArgumentParser(description="Process input folder.")
parser.add_argument("--folder", type=str, help="Path to the input folder (positional)")

args = parser.parse_args()
folder = args.folder


with open(f'{folder}/rasterized_IJ.pkl', 'rb') as f:
    XJ, YJ, J_raw,XI,YI,I, J,muA,muB = pickle.load(f)

I = I.astype(np.float64)

# Add blank padding to the H&E image
padding_width = 300
I_padded = np.pad(
    I,
    pad_width=((0, 0), (padding_width, padding_width), (padding_width, padding_width)),
    mode='constant',
    constant_values=1
)
print(f"original shape of I is : {I.shape}")
print(f"After padding, the shape of I_padded is: {I_padded.shape}")
YI = np.array(range(I_padded.shape[1]))*1. # needs to be longs not doubles for STalign.transform later so multiply by 1.
XI = np.array(range(I_padded.shape[2]))*1. # needs to be longs not doubles for STalign.transform later so multiply by 1.


extentI = STalign.extent_from_x((YI,XI))
extentJ = STalign.extent_from_x((YJ,XJ))
# read from file
pointsJlist = np.load(f'{folder}/SP_points_points.npy', allow_pickle=True).tolist()
print(pointsJlist)
pointsIlist = np.load(f'{folder}/tissue_hires_image_points.npy', allow_pickle=True).tolist()
print(pointsIlist)

pointsIlist_padded = {'p': []}

# Iterate over each point in the original dictionary and add the offset
points_key = pointsIlist.keys()
points_key = list(points_key)[0]
for point in pointsIlist[points_key]:
    x, y = point
    adjusted_x = x + padding_width
    adjusted_y = y + padding_width
    pointsIlist_padded['p'].append((adjusted_x, adjusted_y))


# convert to array
pointsI = []
pointsJ = []


for i in pointsIlist_padded.keys(): # Extract coordinates from pointsIlist_padded
    for j in range(len(pointsIlist_padded[i])):
        pointsI.append([pointsIlist_padded[i][j][1], pointsIlist_padded[i][j][0]])
for i in pointsJlist.keys():
    for j in range(len(pointsJlist[i])):
        pointsJ.append([pointsJlist[i][j][1], pointsJlist[i][j][0]])

pointsI = np.array(pointsI)
pointsJ = np.array(pointsJ)

# now arrays
print(pointsI)
print(pointsJ)

# plot
fig,ax = plt.subplots(1,2)
ax[0].imshow((I_padded.transpose(1,2,0).squeeze()), extent=extentI)
ax[1].imshow((J.transpose(1,2,0).squeeze()), extent=extentJ)

trans_offset_0 = mtransforms.offset_copy(ax[0].transData, fig=fig,
                                       x=0.05, y=-0.05, units='inches')
trans_offset_1 = mtransforms.offset_copy(ax[1].transData, fig=fig,
                                       x=0.05, y=-0.05, units='inches')

ax[0].scatter(pointsI[:,1],pointsI[:,0], c='red', s=10)
ax[1].scatter(pointsJ[:,1],pointsJ[:,0], c='red', s=10)

for i in pointsIlist_padded.keys():
    for j in range(len(pointsIlist_padded[i])):
        ax[0].text(pointsIlist_padded[i][j][0], pointsIlist_padded[i][j][1],f'{i}{j}', c='red', transform=trans_offset_0, fontsize= 8)
for i in pointsJlist.keys():
    for j in range(len(pointsJlist[i])):
        ax[1].text(pointsJlist[i][j][0], pointsJlist[i][j][1],f'{i}{j}', c='red', transform=trans_offset_1, fontsize= 8)

ax[0].set_title('source with pointsI after padding', fontsize=15)
ax[1].set_title('target with pointsJ', fontsize=15)

# invert only rasterized image
ax[0].invert_yaxis()
plt.savefig(f'{folder}/1.2.I_&_S_points.png')

# set device for building tensors
if torch.cuda.is_available():
    torch.set_default_device('cuda:0')
else:
    torch.set_default_device('cpu')

torch.__version__

# compute initial affine transformation from points
L,T = STalign.L_T_from_points(pointsI,pointsJ)
A = STalign.to_A(torch.tensor(L),torch.tensor(T))

print(I_padded.dtype)
print(J.dtype)
print(XI.dtype)
print(YI.dtype)
print(XJ.dtype)
print(YJ.dtype)
# Convert I to a PyTorch Float tensor
print('XI:',type(XI))
print('A:',type(A))
print('I:',type(I_padded))
print('A:',A)
print('YJ:',YJ)
print('XI:',XI)

# compute initial affine transformation from points
AI= STalign.transform_image_source_with_A(A, [YI,XI],I_padded, [YJ,XJ])

#switch tensor from cuda to cpu for plotting with numpy
if AI.is_cuda:
    AI = AI.cpu()

fig,ax = plt.subplots(1,2)

ax[0].imshow((AI.permute(1,2,0).squeeze()), extent=extentJ)
ax[1].imshow((J.transpose(1,2,0).squeeze()), extent=extentJ)

ax[0].set_title('source with affine transformation', fontsize=15)
ax[1].set_title('target', fontsize=15)
plt.savefig(f'{folder}/1.3.source with affine transformation.png')

# run LDDMM
# specify device (default device for STalign.LDDMM is cpu)
if torch.cuda.is_available():
    device = 'cuda:0'
else:
    device = 'cpu'
# Convert muA and muB to tensors with the same shape as J

params = {'L':L,'T':T,
          'pointsI': pointsI,
          'pointsJ': pointsJ,
          'a': 25.0,
        'niter': 8000,   
        'diffeo_start':5000,   
        'epL': 5e-10,               
        'epT': 5e-3,              
        'epV': 0.3, 
        'sigmaM': 0.3, 
        'sigmaB': 0.1, 
        'sigmaA': 0.3,             
        'sigmaP':0.5,          
        'sigmaR': 8e5,          
        'device': 'cuda' if torch.cuda.is_available() else 'cpu',  #
        'dtype': torch.float32,
        'muB': torch.tensor(muB), # black is background in target
        'muA':torch.tensor(muA),
        }

# Run LDDMM and save the resulting plots as PDF
# Convert J to grayscale

# Run LDDMM with the grayscale version of J
# Set a random seed for reproducibility
torch.manual_seed(42)
np.random.seed(42)
out = STalign.LDDMM([YI, XI],I_padded, [YJ, XJ], J, **params)

# Save all the resulting plots as PDF
for i, fig in enumerate(plt.get_fignums()):
    plt.figure(fig)
    plt.savefig(f'{folder}/1.5.LDDMM_out_{i + 1}.pdf', dpi=400)
    plt.savefig(f'{folder}/1.5.LDDMM_out_{i + 1}.png', dpi=400)  # Optional: save for detailed inspection

with open(f'{folder}/2.1.out.pkl', 'wb') as f:
    pickle.dump(out, f)
