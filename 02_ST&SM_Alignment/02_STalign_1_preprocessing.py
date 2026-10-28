import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
from STalign import STalign
import pickle
import argparse
import matplotlib.pyplot as plt
from utils import interactive_roi_hist
from utils import boost_contrast_highlights_and_shadows
from utils import display_image_and_save_gray_value 
from utils import match_histograms_with_nan
from skimage.filters import gaussian
# Parse command-line arguments
parser = argparse.ArgumentParser(description="Process input folder.")
parser.add_argument("--folder", type=str, help="Path to the input folder (positional)")
args = parser.parse_args()
folder = args.folder

#df=pd.read_csv(f'{folder}/pixel_tsne.csv')
hs = pd.read_csv(f'{folder}/hotmodule_scores.csv')
select_modules = int(input("Enter the index for the selected module (e.g., 1): "))

# make plots bigger
plt.rcParams["figure.figsize"] = (15,15)                                                                                                                                                                                                
# get cell centroid coordinates
sh=np.array(hs.iloc[:,select_modules])

if os.path.exists(f'{folder}/Jh_temp.pkl'):
    with open(f'{folder}/Jh_temp.pkl', 'rb') as f:
        XJ, YJ, Jh = pickle.load(f)
else:
    XJ, YJ, Jh, fig = STalign.rasterize(xJ, yJ, sh, dx=1)
    with open(f'{folder}/Jh_temp.pkl', 'wb') as f:
        pickle.dump((XJ, YJ, Jh), f)

# get extent of images
J = np.vstack((Jh, Jh, Jh))
extentJ = STalign.extent_from_x((YJ,XJ))
fig, ax = plt.subplots()
ax.imshow(J.transpose(1, 2, 0))
plt.show()
 
image_file=f'{folder}/tissue_hires_image.png'
V = plt.imread(image_file)

# plot
fig,ax = plt.subplots()
ax.imshow(V)
plt.savefig(f'{folder}/0.2.image_show.png')

Vnorm = STalign.normalize(V)
print("The range of V after normalization is {} to {}".format(Vnorm.min(), Vnorm.max() ))

I = Vnorm.transpose(2,0,1)
print("The new shape of I is {}".format(I.shape))
YI = np.array(range(I.shape[1]))*1. # needs to be longs not doubles for STalign.transform later so multiply by 1.
XI = np.array(range(I.shape[2]))*1. # needs to be longs not doubles for STalign.transform later so multiply by 1.
extentI = STalign.extent_from_x((YI,XI))

fig,ax = plt.subplots()
ax.imshow(I.transpose(1, 2, 0), extent=extentI)
plt.show()

J_inverted = 1 - J
# Plot the inverted image
fig, ax = plt.subplots()
ax.imshow(J_inverted.transpose(1, 2, 0), extent=extentJ)
plt.show()

I_binned= np.mean(I, axis=0)

block_size = 4
h, w = I_binned.shape
h_crop = (h // block_size) * block_size
w_crop = (w // block_size) * block_size
I_binned = I_binned[:h_crop, :w_crop]

output_shape = (I_binned.shape[0] // block_size, I_binned.shape[1] // block_size)
I_binned = I_binned.reshape(output_shape[0], block_size, output_shape[1], block_size).mean(axis=(1, 3))

fig, ax = plt.subplots(1, 2, figsize=(12, 6))
ax[0].imshow(I_binned, cmap='gray')
ax[0].set_title('Original Image I')
ax[1].imshow(I_binned, cmap='gray')
ax[1].set_title('Binned Image I (2x2 blocks)')
plt.tight_layout()
plt.savefig(f'{folder}/0.1.blocks_image.png')
plt.show()

I_gray=np.stack([np.mean(I, axis=0), np.mean(I, axis=0), np.mean(I, axis=0)], axis=0)
print("The shape of I_gray is {}".format(I_gray.shape))

I_adjusted_G = gaussian(I_gray, sigma=2)

vals, patch = interactive_roi_hist(I_binned)

# remove spots with intensity over 0.8
J_thresholded = np.copy(J_inverted)
# Plot the original and modified images
fig, ax = plt.subplots(figsize=(12, 6))
ax.hist(J_thresholded.ravel(), bins=256, histtype='step', color='blue', label='J_thresholded')
ax.set_title('Histogram before Thresholding')
plt.show()

# Prompt the user to input a threshold value
threshold_value = float(input("Enter a threshold value (e.g., 0.6): "))
J_thresholded[J_thresholded > threshold_value] = np.nan

fig, ax = plt.subplots(1, 3, figsize=(12, 6))
ax[0].hist(J_inverted.ravel(), bins=256, histtype='step', color='blue', label='J_inverted_enhanced')
ax[0].set_title('Histogram before Thresholding')
ax[1].hist(J_thresholded.ravel(), bins=256, histtype='step', color='blue', label='J_thresholded')
ax[1].set_title('Histogram after Thresholding')
ax[2].imshow(J_thresholded.transpose(1,2,0).squeeze())
ax[2].set_title('Image after Thresholding (Threshold=0.8)')
ax[2].axis('off')
plt.tight_layout()
plt.show()


PATCH= np.stack([patch, patch, patch])

# Use skimage.exposure.match_histograms to match the pixel distribution of J_thresholded (the spatial metabolomics image) to the H&E image
J_enhanced = match_histograms_with_nan(J_thresholded,PATCH)
J_invert=(1-J_enhanced)
# Plot the adjusted histogram
fig,ax = plt.subplots(1,3)
# Original image histograms
ax[0].hist(J_thresholded.ravel(), bins=256, histtype='step', color='blue', label='J_thresholded')
ax[0].hist(PATCH.ravel(), bins=256, histtype='step', color='red', label='I_gray')
ax[0].set_title('Histograms Before Matching')
ax[0].legend()
# After matching histograms
ax[1].hist(PATCH.ravel(), bins=256, histtype='step', color='red', label='I_gray')
ax[1].hist(J_enhanced.ravel(), bins=256, histtype='step', color='green', label='J_enhanced')
ax[1].set_title('Histogram After Matching')
ax[2].imshow(J_enhanced.transpose(1,2,0).squeeze(), cmap='gray')
plt.tight_layout()
plt.savefig(f'{folder}/0.2.hist_match.png')
plt.show()

#saved_gray_value = display_image_and_save_gray_value(I_adjusted_G,PATCH, J_enhanced, J_invert)

nan_mask = np.isnan(J_enhanced)
# Replace NaN regions with noise
J_enhanced[nan_mask] = 1

J_enhanced_inverted = 1 - J_enhanced # Invert again so that the sample is bright and the background is dark
fig, ax = plt.subplots(figsize=(12, 6))
ax.imshow(J_enhanced_inverted.transpose(1,2,0).squeeze())
plt.tight_layout()
plt.show()

# Apply piecewise linear enhancement to shadows and highlights, leaving midtones unchanged
# For input pixels below low_thresh, set the value to gain_low * x and clip to (0, 1)
# For input pixels above high_thresh, set the value to high_thresh + gain_high * (x - high_thresh) and clip to (0, 1)
# Set low_thresh to 0 below to enhance only bright regions, leaving shadows unchanged
J_enhanced_inverted_boost = boost_contrast_highlights_and_shadows(
    J_enhanced_inverted, low_thresh=0, high_thresh=0.3, gain_low=0, gain_high=5
)

# Display the boosted image
fig, ax = plt.subplots(figsize=(12, 6))
ax.imshow(J_enhanced_inverted_boost.transpose(1, 2, 0).squeeze(), cmap='gray')
ax.set_title('Boosted Image')
ax.axis('off')
plt.tight_layout()
plt.savefig(f'{folder}/0.3.boosted_image.png')
plt.show()


muA=display_image_and_save_gray_value(I_adjusted_G, PATCH, J_enhanced, J_enhanced_inverted_boost)
muA=[muA,muA,muA]
print(muA)
muB=[0,0,0]
print(muB)

I = I.astype(np.float64)  # Convert I to float64
with open(f'{folder}/rasterized_IJ.pkl', 'wb') as f:
    pickle.dump((XJ,YJ,J,XI,YI,I,J_enhanced_inverted_boost,muA,muB), f)

from os.path import basename

import os
np.savez(f'{folder}/SP_points', x=XJ,y=YJ,I=J)
np.savez(f'{folder}/tissue_hires_image', x=XI,y=YI,I=I)    


# savelog
log_file = f'{folder}/processing_log.txt'
with open(log_file, 'w') as f:
    f.write("Processing Log\n")
    f.write("=====================\n")
    f.write(f"Input folder: {folder}\n")
    f.write(f"Selected module index (for hot module): {select_modules}\n")
    f.write(f"Threshold value for J_thresholded: {threshold_value}\n")
    f.write(f"Gaussian smoothing sigma for I_adjusted_G: 2\n")
    f.write(f"Boost contrast settings:\n")
    f.write(f"  - low_thresh: 0\n")
    f.write(f"  - high_thresh: 0.3\n")
    f.write(f"  - gain_low: 0\n")
    f.write(f"  - gain_high: 5\n")
    f.write("\nGenerated Files:\n")
    f.write(f"  - {folder}/0.2.image_show.png\n")
    f.write(f"  - {folder}/rasterized_IJ.pkl\n")
    f.write(f"  - {folder}/J_tsne.pkl\n")
    f.write(f"  - {folder}/{basename(folder)}.npz\n")
    f.write(f"  - {folder}/tissue_hires_image.npz\n")
    f.write("=====================\n")
