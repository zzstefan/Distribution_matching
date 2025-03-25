% Transforms the computed tracts into the voxel space.
%
% [X, seedPointsVox, ind, mask] = tracts2vox(tracts, mask)
%  
% tracts:         Tract structure resulting from HoughTract.m.
%                 (It can also be a cell array; see 'X' below.)
% mask:           Optional mask(s), either a single 3D image or a cell of
%                 multiple 3D images, each containing a region of interest.
%                 If provided, results will be limited to only the tracts
%                 that intersect all the masks.
%
% X:              A cell array of the selected tracts, with each element
%                 being a tall matrix of the coordinates of the tract
%                 points in the image voxel space.
%                 (If 'X' is fed back as the input parameter 'tracts', the
%                 function will use it to only run the masking part.)
% seedPointsVox:  Coordinates of the corresponding seed points, in the
%                 image voxel space.
% ind:            Indices of the chosen tracts (that pass through all the
%                 masks). For example, the score of the tract X{n} is
%                 tracts.score(ind(n)). 'ind' can be used with
%                 sampleTracts.m to make a new tract structure with only
%                 the chosen tracts.
%
% See also:  HoughTract, genSeedPoints, showTracts, exportTrackVis,
%            exportText, sampleTracts, connMatrix, fiberDensityVol,
%            CSAODF_CLI, EXAMPLE, EXAMPLE_CLI.

% Code by Iman Aganj.

function [X, seedPointsVox, ind, mask] = tracts2vox(tracts, mask)

if isstruct(tracts)
    [seedPoints, polyCoeffs, partialLengths, box2vox, order, resS, maxS] = deal(tracts.seedPoints, tracts.polyCoeffs, tracts.partialLengths, tracts.box2vox, tracts.param.order, tracts.param.resS, tracts.param.maxS);
    seqS = 0:resS:maxS;
    nSeqS = length(seqS);
    nSeedPoints = size(seedPoints, 1);
    xCurve = zeros(nSeqS, 3, 2);
    for iDir = 1:2
        sPower(:,:,iDir) = bsxfun(@power, (2*iDir-3)*seqS, (0:order-1)');
    end
    X = cell(nSeedPoints, 1);
    for i = 1:nSeedPoints
        for iDir = 1:2
            theta = polyCoeffs(i,1:order) * sPower(:, :, iDir);
            phi = polyCoeffs(i,order+1:2*order) * sPower(:, :, iDir);
            sTheta = sin(theta);
            cTheta = cos(theta);
            xCurve(:,:,iDir) = [seedPoints(i,:); ((2*iDir-3)*resS)*[sTheta(1:end-1).*cos(phi(1:end-1)); sTheta(1:end-1).*sin(phi(1:end-1)); cTheta(1:end-1)]'];
        end
        x = cumsum(xCurve);
        x = [x(partialLengths(i,1):-1:2,:,1); x(1:partialLengths(i,2),:,2)];
        X{i} = [x ones(size(x,1),1)] * box2vox';
    end
    seedPointsVox = [seedPoints ones(nSeedPoints,1)] * box2vox';
elseif iscell(tracts)
    X = tracts;
    nSeedPoints = length(X);
    seedPointsVox = nan(nSeedPoints, 3);
else
    error('''tracts'' must be either a structure or a cell array.')
end

if exist('mask', 'var') && ~isempty(mask)
    if ~iscell(mask)
        mask = mat2cell(mask, size(mask,1), size(mask,2), size(mask,3));
    end
    sz = size(mask{1});
    ind = true(nSeedPoints, 1);
    for n = 1:nSeedPoints
        x = int32(X{n});
        legitX = all(x>0, 2) & all(bsxfun(@le, x, sz),2);
        for i = 1:length(mask)
            ind(n) = ind(n) && any(mask{i}(sub2ind(sz, x(legitX,1), x(legitX,2), x(legitX,3))));
        end
    end
    X = X(ind);
    seedPointsVox = seedPointsVox(ind,:);
    ind = find(ind);
end
