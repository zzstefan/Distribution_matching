% Exports the computed tracts to a .trk file, so it can be visualized by
% TrackVis (www.trackvis.org) and Freeview (https://surfer.nmr.mgh.harvard.edu).
% An optional parameter will be ignored if [] is entered for it. 
%
% exportTrackVis(filename, tracts, nTracts, mask, vox2ras, voxelSize, voxelOrder, imgOrient)
%  
% filename:   Output '.trk' file name, where tracts (sorted, high-score
%             first) are exported.
% tracts:     Tract structure resulting from HoughTract.m.
% nTracts:    Optional number of top-scored tracts to be exported.
% mask:       Optional mask(s), either a single 3D image or a cell of
%             multiple 3D images, each containing a region of interest. If
%             provided, only the tracts that intersect all of the masks will
%             be exported.
% vox2ras:    Optional 3x4 (or 4x4) vox2ras matrix.
% voxelSize:  Optional 3x1 voxel size.
% voxelOrder: Optional voxel order (default: 'LPS').
% imgOrient:  Optional 3x2 image orientation (default: [1 0; 0 1; 0 0]).
%
% See also:   HoughTract, exportText, tracts2vox, showTracts, sampleTracts,
%             connMatrix, fiberDensityVol, CSAODF_CLI, EXAMPLE, EXAMPLE_CLI.

% Code by Iman Aganj.

function exportTrackVis(filename, tracts, nTracts, mask, vox2ras, voxelSize, voxelOrder, imgOrient)

score = tracts.score;
if exist('mask', 'var') && ~isempty(mask)
    [X, ~, ind] = tracts2vox(tracts, mask);
    score = score(ind);
else
    X = tracts2vox(tracts);
end
if ~exist('nTracts', 'var') || isempty(nTracts)
    nTracts = length(X);
end
if exist('vox2ras', 'var') && ~isempty(vox2ras)
    if all(size(vox2ras)==[3 4])
        vox2ras = [vox2ras; 0 0 0 1];
    elseif all(size(vox2ras)==[4 4])
        if any(vox2ras(4,:)~=[0 0 0 1])
            warning('vox2ras(4,:) is expected to be [0 0 0 1]. Did you mean to input vox2ras''?')
        end
    else
        error('''vox2ras'' must be either 4x4 or 3x4!')
    end
else
    vox2ras = zeros(4);
end
if ~exist('voxelOrder', 'var') || isempty(voxelOrder)
    voxelOrder = 'LPS';
end
if exist('voxelSize', 'var') && ~isempty(voxelSize)
    if length(voxelSize) ~= 3
        error('voxelSize must have length of 3.')
    end
else
    voxelSize = diag(tracts.box2vox); % Needs modification
    voxelSize = max(voxelSize)./ voxelSize;
end
if ~exist('imgOrient', 'var') || isempty(imgOrient)
    imgOrient = [1 0; 0 1; 0 0]; % image_orientation_patient?
end

[score, ind] = sort(score, 'descend');
X = X(ind);
fid = fopen(filename, 'wb');

% Header
fwrite(fid, ['TRACK' 0], 'char');
fwrite(fid, tracts.param.sz(1:3), 'short');
fwrite(fid, voxelSize, 'float');
fwrite(fid, [0 0 0], 'float');
fwrite(fid, 1, 'short');
propertyName = 'Score';
fwrite(fid, [propertyName blanks(200-length(propertyName))], 'char');
fwrite(fid, 0, 'short');
fwrite(fid, blanks(200), 'char');
fwrite(fid, vox2ras', 'float');
fwrite(fid, blanks(444), 'char');
fwrite(fid, [voxelOrder 0], 'char'); % voxel_order needs verification
fwrite(fid, [voxelOrder 0], 'char'); % pad2?
fwrite(fid, imgOrient, 'float');
fwrite(fid, [0 0], 'char'); % pad1?
fwrite(fid, 0, 'uchar');
fwrite(fid, 0, 'uchar');
fwrite(fid, 0, 'uchar');
fwrite(fid, 0, 'uchar');
fwrite(fid, 0, 'uchar');
fwrite(fid, 0, 'uchar');
fwrite(fid, nTracts, 'int');
fwrite(fid, 2, 'int');
fwrite(fid, 1000, 'int');

% Tracts
for n = 1:nTracts
    nPts = size(X{n},1);
    fwrite(fid, nPts, 'int');
%     X{n} = [X{n} ones(nPts,1)] * vox2ras';
%     X{n} = X{n}(:,1:3);
    fwrite(fid, [(voxelSize(:)'.*X{n})'; score(n)*ones(1, nPts)], 'float');
end

fclose(fid);
