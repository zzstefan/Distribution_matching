% This is an example script showing how to use this toolbox to compute
% q-ball CSA-ODFs from diffusion-weighted MRI, perform Hough-transform
% tractography, and visualize the ODFs, tracts, and connectivity matrix.
% For a user-friendly Matlab-compilable pipeline, see EXAMPLE_CLI.m.
%
% CSA-ODF is based on the method proposed in:
%
% I. Aganj, C. Lenglet, G. Sapiro, E. Yacoub, K. Ugurbil, and N. Harel,
% "Reconstruction of the orientation distribution function in single and
% multiple shell q-ball imaging within constant solid angle,"
% Magnetic Resonance in Medicine, vol. 64, no. 2, pp. 554-566, 2010.
%
% The Hough-transform approach to fiber tractography has been proposed in:
%
% I. Aganj, C. Lenglet, N. Jahanshad, E. Yacoub, N. Harel, P. Thompson, and
% G. Sapiro, "A Hough transform global probabilistic approach to
% multiple-subject diffusion MRI tractography,"
% Medical Image Analysis, vol. 15, no. 4, pp. 414-425, 2011.
%
% The gradient-table verification approach was introduced in:
%
% I. Aganj, "Automatic verification of the gradient table in
% diffusion-weighted MRI based on fiber continuity,"
% Scientific Reports, vol. 8, Article no. 16541, 2018.
%
% The connectivity matrix augmentation approach was introduced in:
%
% I. Aganj, G. Prasad, P. Srinivasan, A. Yendiki, P. M. Thompson, and B.
% Fischl, "Structural brain network augmentation via Kirchhoff's laws,"
% in Proc. Joint Annual Meeting of ISMRM-ESMRMB, Milan, Italy, 2014.
%
% See also:   EXAMPLE_CLI.m

% Codes by Iman Aganj
% http://iman.mgh.harvard.edu

%% Step 1: Download the data
% In this example, we use a subject from the publicly available FreeSurfer
% diffusion tutorial dataset, which is downloaded and extracted here.
% Alternatively, to process your own data, skip this and go to Step 3:
gunzip https://surfer.nmr.mgh.harvard.edu/pub/data/subjects/elmo.2012.tar.gz
untar elmo.2012.tar
delete elmo.2012.tar

%% Step 2: Read the data
% The "tutorial_data" folder is assumed to exist in the current folder.
% To process your own data instead, skip this and go to Step 3.
subName = 'elmo.2012';

% Read the gradient table.
bvals = load(fullfile(subName,'dmri','bvals'), 'ascii');
bvecs = load(fullfile(subName,'dmri','bvecs'), 'ascii');
b = unique(bvals);
disp(array2table([b sum(bsxfun(@eq, b(:,1), bvals'),2)], 'VariableNames', {'bValue', 'numberOfImages'}))
bval = b(2);
indB0 = bvals==0;
indB  = bvals==bval;
bvecs = bvecs(indB,:);

% Read the images
dataFile = fullfile(subName,'dmri','dwi.nii');
if isunix, dataFile = [dataFile '.gz']; else, gunzip([dataFile '.gz']), end
V = load_nifti(dataFile); % load_nifti.m is part of FreeSurfer's Matlab codes (https://surfer.nmr.mgh.harvard.edu). Matlab's 'niftiread' can also be used.
vox2ras = [V.srow_x V.srow_y V.srow_z]';
voxelSize = V.pixdim(2:4);
S = V.vol(:,:,:,indB); % diffusion-weighted images
S0 = mean(V.vol(:,:,:,indB0),4); % b=0 image

% Create a mask
segFile = fullfile(subName,'dlabel','diff','aparc+aseg.bbr.nii');
if isunix, segFile = [segFile '.gz']; else, gunzip([segFile '.gz']), end
V = load_nifti(segFile);
meanADC = (log(S0) - mean(log(S),4)) / bval;
Mask = V.vol>0 & meanADC<1e-3;
V.vol = Mask;
save_nifti(V, fullfile(subName,'Mask.nii'));
clear V

%% Step 3: Verify the gradient table
% It is assumed that S contains the diffusion images as a 4D array (or a 1D
% cell array of 3D images), S0 is the non-diffusion-weighted image, bvecs
% is the gradient table (tall matrix), and Mask is the brain mask excluding
% the CSF.

% Make sure the gradient table is not misaligned with the image due to a
% permutation and/or a flip of the coordinates. If it is, correct the
% b-vectors.
[~, ~, correctedBVecs] = checkGradTable(S, S0, bvecs, Mask, 0, voxelSize);

%% Step 4: Reconstruct the CSA-ODFs
% Compute the gradient table in spherical coordinates.
angles = cart2sph_phys(bvecs);

% Compute the ODFs. Parameters such as the spherical harmonic basis order,
% delta (regularization), and indexGPU (for faster computation) can be
% optionally specified (see the help of reconCSAODF.m). To reconstruct the
% ODFs from three shells (with the same gradient directions for all shells,
% and b-values that are an arithmetic sequence), see reconCSAODF3Q.m.
sh = reconCSAODF(S, S0, angles);

% Compute the Generalized Fractional Anisotropy.
GFA = makeGFA(sh);

%% Step 5: Visualize the ODFs
% Visualize some of the ODFs to verify the results. ODFs can also be
% sampled in a set of directions using sampleODFs.m.
ROI = round(bsxfun(@times, size(GFA)', [.3 .7; .5 .5; .2 .9]));
showODFs(sh, ROI, GFA)

%% Step 6: Hough-transform tractography
clear param

% Approximate the scalar fiber probability map as the masked GFA.
pVol = GFA .* Mask;

% Create the transformation from the tractography box to the image. If
% vox2ras isn't available, a 1x3 vector of pixel sizes will also work.
box2vox = makeBox2vox(size(pVol), vox2ras);

% Generate 10000 seed points.
seedPoints = genSeedPoints(10000, pVol, box2vox);

% Choose the suggested tract length prior. Increasing/decreasing lambda
% results in longer/shorter tracts, respectively.
param.lambda = suggestLambda(sh, pVol);

% Ask for the maximum number of available cores to be used. If a GPU is
% available, using the param.indexGPU parameter instead can significantly
% increase speed. There are additional parameters that change the
% precision-speed tradeoff (see the help of HoughTract.m).
param.nCores = inf;

% Perform Hough-transform tractography. This step can take hours depending
% on the chosen parameter values. For a quick test, choose a low param.order=2.
tracts = HoughTract(sh, pVol, seedPoints, box2vox, param);
save(fullfile(subName, 'tracts'), 'tracts')

%% Step 7: Visualize the tracts
% Visualize the tracts of the 6000 top-scored seed points on the GFA volume.
showTracts(tracts, GFA, 6000)

% Export the 6000 top-scored tracts to a .trk file, so it can be visualized
% by TrackVis (www.trackvis.org) and Freeview (https://surfer.nmr.mgh.harvard.edu).
% Similarly, exportText.m exports the tract coordinates into a text file. 
exportTrackVis(fullfile(subName, 'Tracts.trk'), tracts, 6000, [], vox2ras, voxelSize);

% Visualize the tracts passing through the thalamus, along with their seed
% points. The FreeSurfer color look-up table is in FreeSurferColorLUT.txt.
segVol = load_nifti(segFile); segVol = segVol.vol;
ROI = segVol==10 | segVol==49;
showTracts(tracts, GFA, [], ROI, [], true)

% Get the coordinates (X) and indices (ind) of the tract points going
% through both the left thalamus (ROI1) and the left-hemisphere cortex (ROI2).
ROI1 = segVol==10;
ROI2 = segVol>=1000 & segVol<2000;
[X, ~, ind] = tracts2vox(tracts, {ROI1, ROI2});

% Make a fiber density volume from the chosen tracts and plot its isosurface.
densVol = fiberDensityVol(sampleTracts(tracts, ind));
figure, isosurface(densVol, 0), axis equal tight

%% Step 8: Compute and visualize the connectivity matrix
% Choose labels.
labels = [8 10:13 17:18 28 1001:1003 1005:1035 16 47 49:54 60 2001:2003 2005:2035]';
labelNames = FreeSurferLUT(labels);

% Create the connectivity matrix by counting the number (weighted by score)
% of tracts passing through every pair of ROIs. Then visualize the matrix
% and the tracts interactively; moving the mouse pointer on the matrix
% updates the tractogram showing only the tracts connecting the
% corresponding pair of ROIs. (Clicking will lock the pointer's position.)
C = connMatrix(tracts, segVol, labels, true, true, labelNames);

% To obtain a fuller matrix, use a larger param.lambda in tractography.

%% Step 9: Augment the connectivity matrix
% Augment the connectivity matrix with indirect connections by modeling it
% as electric conductance.
aC = augConnMatrix(C);

% Show the matrices side by side.
figure
for i = 1:2
    subplot(1,2,i)
    if i==1
        imagesc(C)
        title('Original connectivity matrix')
    else
        imagesc(aC)
        title('Augmented connectivity matrix')
    end
    axis equal tight
    colormap(jet(2^12))
    colorbar
    xticks(1:size(C,2)), xticklabels(labelNames)
    yticks(1:size(C,1)), yticklabels(labelNames)
end
