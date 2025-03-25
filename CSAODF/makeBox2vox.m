% Computes a 3x4 affine transformation matrix that takes the seed points
% from the [-1,1]^3 box to the image voxel space.
%
% box2vox = makeBox2vox(sz, pixdimORvox2ras)
%  
% sz:               1x3 vector of image size.
% pixdimORvox2ras:  Either a 3x4 (or 3x3) vox2ras matrix, or a 1x3 vector
%                   of voxel dimensions.
%
% box2vox:          Computed transformation.
%
% See also:  HoughTract, genSeedPoints, tracts2vox, exportTrackVis,
%            exportText, sampleTracts, fiberDensityVol, CSAODF_CLI,
%            EXAMPLE, EXAMPLE_CLI.

% Code by Iman Aganj.

function box2vox = makeBox2vox(sz, pixdimORvox2ras)

if ~(all(size(pixdimORvox2ras)==[3 4]) || all(size(pixdimORvox2ras)==[3 3]) || all(size(pixdimORvox2ras)==[1 3]) || all(size(pixdimORvox2ras)==[3 1]))
    error('pixdimORvox2ras must be either a 3x4 (or 3x3) vox2ras matrix, or a 3x1 vector of voxel dimensions.')
end
if ~all(size(sz)==[1 3])
    error('sz must be the 1x3 vector of image size.')
end
if min(size(pixdimORvox2ras))>1
    pixdim = svd(pixdimORvox2ras(1:3,1:3))';
else
    pixdim = pixdimORvox2ras(:)';
end

box2vox = [diag(max(((sz-1)/2).*pixdim) ./ pixdim) (sz'+1)/2];
