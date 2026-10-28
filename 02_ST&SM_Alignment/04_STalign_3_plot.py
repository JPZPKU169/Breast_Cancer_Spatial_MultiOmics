import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import torch
from STalign import STalign
import pickle
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import argparse

parser = argparse.ArgumentParser(description="Process input folder.")
parser.add_argument("--folder", type=str, help="Path to the input folder (positional)")

args = parser.parse_args()
folder = args.folder

# Load the image
with open(f'{folder}/rasterized_IJ.pkl', 'rb') as f:
    XJ, YJ, J_raw,XI,YI,I, J,muA,muB = pickle.load(f) 

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

# Load STalign results
with open(f'{folder}/2.1.out.pkl', 'rb') as f:
    out=pickle.load(f)

# Load landmarks
pointsJlist = np.load(f'{folder}/SP_points_points.npy', allow_pickle=True).tolist()
pointsIlist = np.load(f'{folder}/tissue_hires_image_points.npy', allow_pickle=True).tolist()

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

L_pre,T_pre = STalign.L_T_from_points(pointsI,pointsJ)
A_pre = STalign.to_A(torch.tensor(L_pre),torch.tensor(T_pre))

df1 = pd.read_csv(f'{folder}/coords.csv')
print(df1.head())
xI = np.array(df1['x'])
yI = np.array(df1['y'])

#### get necessary output variables
A_post = out['A']
v = out['v']
xv = out['xv']
WM = out['WM'] 
WB=out['WB']
WA=out['WA']
print('A_pre:',A_pre)

print('A_post:',A_post)
if torch.cuda.is_available():
    device = 'cuda:0'
else:
    device = 'cpu'

# Convert all tensors to float32 and move them to the same device
YI = torch.tensor(YI, dtype=torch.float32, device=device)
XI = torch.tensor(XI, dtype=torch.float32, device=device)
I_padded = torch.tensor(I_padded, dtype=torch.float32, device=device)
YJ = torch.tensor(YJ, dtype=torch.float32, device=device)
XJ = torch.tensor(XJ, dtype=torch.float32, device=device)
J = torch.tensor(J, dtype=torch.float32, device=device)

A_pre = torch.tensor(A_pre, dtype=torch.float32, device=device)
A_post = torch.tensor(A_post, dtype=torch.float32, device=device)

# Compute the source image before and after the linear transformation
AI_pre= STalign.transform_image_source_with_A(A_pre, [YI,XI], I_padded, [YJ,XJ]) # Linear transformation before LDDMM 
AI_post= STalign.transform_image_source_with_A(A_post, [YI,XI], I_padded, [YJ,XJ]) # LDDMM linear transformation

# Compute the affine transformation of all corresponding points
affine_spot = np.matmul(np.array(A_post.cpu()),np.array([yI, xI, np.ones(len(xI))]))
xIaffine = affine_spot[1,:]
yIaffine = affine_spot[0,:]
# Compute the transformed positions of all landmarks
ypointsI = pointsI[:,0]
xpointsI = pointsI[:,1]
affine_pre = np.matmul(np.array(A_pre.cpu()),np.array([ypointsI, xpointsI, np.ones(len(ypointsI))])) 
affine_post = np.matmul(np.array(A_post.cpu()),np.array([ypointsI, xpointsI, np.ones(len(ypointsI))]))
pointsIaffine_pre = np.column_stack((affine_pre[0,:],affine_pre[1,:]))
pointsIaffine_post = np.column_stack((affine_post[0,:],affine_post[1,:]))

# Apply linear and nonlinear transformations to the source
pointsI=torch.tensor(pointsI, dtype=torch.float32, device=device)
# apply transform
phii = STalign.build_transform(xv,v,A_post,XJ=[YJ,XJ],direction='b') # Compute the deformation grid for the combined linear and nonlinear transformation
phiI = STalign.transform_image_source_to_target(xv,v,A_post,[YI,XI],I_padded,[YJ,XJ]) # Apply the nonlinear transformation to the source image
phiipointsI = STalign.transform_points_source_to_target(xv,v,A_post,pointsI) # Apply the nonlinear transformation to landmarks

with open(f'{folder}/3.1.phiI.pkl', 'wb') as f:
    pickle.dump((extentI,phii,phiI),f)

#switch tensor from cuda to cpu for plotting with numpy
if phii.is_cuda:
    phii = phii.cpu()
if phiI.is_cuda:
    phiI = phiI.cpu()
if XJ.is_cuda:
    XJ = XJ.cpu()
if YJ.is_cuda:
    YJ = YJ.cpu()
if phiipointsI.is_cuda:
    phiipointsI = phiipointsI.cpu()
if J.is_cuda:
    J = J.cpu()

# Deformation grid for the nonlinear transformation
fig,ax = plt.subplots(figsize=(20,20))
levels = np.arange(-2000,5000,25) # The last argument is the step size
ax.contour(XJ,YJ,phii[...,0],colors='r',linestyles='-',levels=levels)
ax.contour(XJ,YJ,phii[...,1],colors='g',linestyles='-',levels=levels)
ax.set_aspect('equal')
ax.set_title('source to target')
ax.imshow(phiI.permute(1,2,0)/torch.max(phiI),extent=extentJ)
ax.scatter(phiipointsI[:,1].detach(),phiipointsI[:,0].detach(),c="m")
fig.savefig(f'{folder}/2.1.source_elastic.png',dpi=900)

# Before LDDMM
fig,ax = plt.subplots()
ax.imshow(J.permute(1, 2, 0), extent=extentJ)
I_np = AI_pre.detach().cpu().numpy() if hasattr(AI_pre, 'detach') else AI_pre
# Transpose dimensions from (C, H, W) to (H, W, C)
I_np = I_np.transpose(1, 2, 0).squeeze()
# Display the image and set transparency
ax.imshow(I_np, extent=extentJ, alpha=0.4)  # Set alpha to 0.5 to increase transparency
ax.scatter(pointsIaffine_pre[:, 1], pointsIaffine_pre[:, 0], c="blue", label='source landmarks aligned', s=5)
ax.scatter(pointsJ[:, 1], pointsJ[:, 0], c='red', label='target landmarks', s=5)
ax.set_aspect('equal')
lgnd = plt.legend(loc="upper right", scatterpoints=1, fontsize=10)
for handle in lgnd.legend_handles:
    handle.set_sizes([10.0])
ax.set_title('pre_LDDMM', fontsize=15)
plt.savefig(f'{folder}/2.2.overlay_pre_LDDMM.png',dpi=600)

# Before the nonlinear transformation
fig,ax = plt.subplots()
ax.imshow(J.permute(1, 2, 0), extent=extentJ)
I_np = AI_post.detach().cpu().numpy() if hasattr(AI_post, 'detach') else AI_post
# Transpose dimensions from (C, H, W) to (H, W, C)
I_np = I_np.transpose(1, 2, 0).squeeze()
# Display the image and set transparency
ax.imshow(I_np, extent=extentJ, alpha=0.4)  # Set alpha to 0.5 to increase transparency
ax.scatter(pointsIaffine_post[:, 1], pointsIaffine_post[:, 0], c="blue", label='source landmarks aligned', s=5)
ax.scatter(pointsJ[:, 1], pointsJ[:, 0], c='red', label='target landmarks', s=5)
ax.set_aspect('equal')
lgnd = plt.legend(loc="upper right", scatterpoints=1, fontsize=10)
for handle in lgnd.legend_handles:
    handle.set_sizes([10.0])
ax.set_title('pre_nonlinear', fontsize=15)
plt.savefig(f'{folder}/2.3.overlay_pre_nonlinear.png',dpi=600)

# After the nonlinear transformation
fig,ax = plt.subplots()
ax.imshow(J.permute(1, 2, 0), extent=extentJ)
phiI_np = phiI.detach().cpu().numpy() if hasattr(phiI, 'detach') else phiI
phiI_np = phiI_np.transpose(1, 2, 0).squeeze()
ax.imshow(phiI_np, extent=extentJ, alpha=0.4)  # Set alpha to 0.5 to increase transparency
ax.scatter(phiipointsI[:, 1], phiipointsI[:, 0], c="yellow", label='source landmarks aligned', s=5)
ax.scatter(pointsJ[:, 1], pointsJ[:, 0], c='red', label='target landmarks', s=5)
ax.set_aspect('equal')
lgnd = plt.legend(loc="upper right", scatterpoints=1, fontsize=10)
for handle in lgnd.legend_handles:
    handle.set_sizes([10.0])
ax.set_title('source and target aligned with aligned pointsI and pointsJ', fontsize=15)
plt.savefig(f'{folder}/2.4.overlay_post_nonlinear.png',dpi=600)
