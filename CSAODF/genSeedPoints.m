% Generates quasi-random seed points for Hough-transform tractography.
%
% [seedPoints, box2vox] = genSeedPoints(nSeedPoints, pVol, box2vox)
%  
% nSeedPoints:  Number of seed points to be generated.
% pVol:         3D probability map of the existence of fiber bundles. For
%               instance, it can be the GFA (that may be created using
%               makeGFA.m) or FA inside of the white+gray matter mask, and
%               0 outside of it.
% box2vox:      An optional 3x4 affine transformation matrix that takes the
%               seed points from the [-1,1]^3 tractography box to the image
%               voxel space, and can be generated using makeBox2vox.m. If
%               not specified, it will be created assuming that the image
%               has isotropic voxels, and will be given in the output.
%
% seedPoints:   nSeedPoints x 3 matrix of the Cartesian coordinates of the
%               seed points (in the cube [-1,1]^3).
%
% See also:  HoughTract, makeBox2vox, showTracts, CSAODF_CLI, EXAMPLE,
%            EXAMPLE_CLI.

% Code by Iman Aganj.

function [seedPoints, box2vox] = genSeedPoints(nSeedPoints, pVol, box2vox)

sz = size(pVol);
if length(sz)<3
    sz = [sz ones(1,3-length(sz))];
end
X = zeros(nSeedPoints,3);
P = haltonset(3);
nHalton = 1;
for i = 1:nSeedPoints
    x = inf(1,3);
    while(any(x<=0) || any(x>sz) || rand(1)>pVol(round(x(1)), round(x(2)), round(x(3))))
        x = P(nHalton,:).*(sz-1) + 1;
        nHalton = nHalton + 1;
    end
    X(i,:) = x;
end
if ~exist('box2vox', 'var')
    box2vox = makeBox2vox(size(pVol), ones(1,3));
end
seedPoints = ([box2vox; 0 0 0 1]\[X ones(nSeedPoints,1)]')';
seedPoints = seedPoints(:,1:3);
