% This is an example script showing how to use the command-line interface
% of this toolbox to compute q-ball CSA-ODFs from diffusion-weighted MRI,
% perform Hough-transform tractography, and visualize the ODFs and the
% tracts. The following commands can be run directly in Matlab, or in the
% system shell after compiling with 'mcc -m CSAODF_CLI.m'. In this example,
% we use a subject from the publicly available FreeSurfer tutorial dataset
% (https://surfer.nmr.mgh.harvard.edu/pub/data/subjects/elmo.2012.tar.gz),
% and assume that we are in the 'elmo.2012' directory. You can process your
% own data similarly. 
% For a more rigorous and flexible Matlab pipeline, see EXAMPLE.m.
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
% multiple-subject diffusion MRI tractography," Medical Image Analysis,
% vol. 15, no. 4, pp. 414-425, 2011.
%
% The gradient-table verification approach was introduced in:
%
% I. Aganj, "Automatic verification of the gradient table in
% diffusion-weighted MRI based on fiber continuity,"
% Scientific Reports, vol. 8, Article no. 16541, 2018.
%
% See also:   EXAMPLE.m

% Codes by Iman Aganj
% http://iman.mgh.harvard.edu


% We assume the 4D diffusion data to be in dmri/dwi.nii.gz, the b-vectors
% in dmri/bvecs, and the b-values in dmri/bvals. Mask.nii is a GM+WM mask,
% which, for this dataset, is generated in EXAMPLE.m.

% (Optional:) Make sure the gradient table is not misaligned with the image
% due to a permutation and/or a flip of the coordinates. Save the corrected
% b-vectors in dmri/corrected_bvecs.
CSAODF_CLI verify dmri/dwi.nii.gz dmri/bvecs dmri/bvals dmri/corrected_bvecs mask Mask.nii

% Reconstruct the CSA-ODFs and save them into ODFs.nii.
% Run 'CSAODF_CLI reconCSAODF' to see all the available options.
CSAODF_CLI reconCSAODF dmri/dwi.nii.gz dmri/bvecs dmri/bvals ODFs.nii

% Compute the Generalized Fractional Anisotropy.
CSAODF_CLI makeGFA ODFs.nii GFA.nii

% Visualize some of the ODFs to verify the results.
% Run 'CSAODF_CLI showODFs' to see all the available options.
CSAODF_CLI showODFs ODFs.nii

% Generate 10000 seed points. Run 'CSAODF_CLI genSeedPoints' to see all the
% available options.
CSAODF_CLI genSeedPoints SeedPoints.mat Mask.nii pVol GFA.nii

% Perform Hough-transform tractography. This step can take hours depending
% on the chosen parameter values. Run 'CSAODF_CLI HoughTract' to see all
% the available options, including how to use a GPU for speed-up, specify
% the number of cores, or change the precision-speed tradeoff. (For a quick
% test, add the 'order 2' option.)
CSAODF_CLI HoughTract ODFs.nii SeedPoints.mat tracts.mat Mask.nii

% Alternatively, for distributed computing:
%     Generate the commands for 10 batch jobs in jobs.txt, each performing
%     a part of the Hough-transform tractography.
CSAODF_CLI HoughTract ODFs.nii SeedPoints.mat tracts.mat Mask.nii createJobs 10 jobFilename jobs.txt
%     After distributing the jobs in jobs.txt and running them, collect and
%     merge the results:
CSAODF_CLI collectHoughResults tracts.mat 10

% Visualize the tracts of the 6000 top-scored seed points on the GFA
% volume. Run 'CSAODF_CLI showTracts' to see all the available options.
CSAODF_CLI showTracts tracts.mat overlay GFA.nii nTracts 6000

% Export the tracts to a .trk file, so it can be visualized by TrackVis
% (www.trackvis.org) and Freeview (https://surfer.nmr.mgh.harvard.edu).
% Similarly, using 'text' instead of 'trk' exports the tract coordinates
% into a text file. Run 'CSAODF_CLI export' to see the details and
% available options.
CSAODF_CLI export trk tracts.mat tracts.trk
