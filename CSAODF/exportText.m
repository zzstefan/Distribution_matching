% Exports the computed tracts into a text file, which will include first
% the number of tracts, and then for each tract the number of points, the
% tract score, and the list of voxel coordinates of the points. An optional
% parameter will be ignored if [] is entered for it.
%
% exportText(filename, tracts, nTracts, mask)
%  
% filename:  Output text file name, where tracts (sorted, high-score first)
%            are exported.
% tracts:    Tract structure resulting from HoughTract.m.
% nTracts:   Optional number of top-scored tracts to be exported.
% mask:      Optional mask(s), either a single 3D image or a cell of
%            multiple 3D images, each containing a region of interest. If
%            provided, only the tracts that intersect all of the masks will
%            be exported.
%
% See also:  HoughTract, exportTrackVis, tracts2vox, showTracts,
%            sampleTracts, connMatrix, fiberDensityVol, CSAODF_CLI,
%            EXAMPLE, EXAMPLE_CLI.

% Code by Iman Aganj.

function exportText(filename, tracts, nTracts, mask)

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
[score, ind] = sort(score, 'descend');
X = X(ind);

fid = fopen(filename, 'w');
fprintf(fid, '%d\n', nTracts);
for n = 1:nTracts
    fprintf(fid, '%d\n%f\n', size(X{n},1), score(n));
    fprintf(fid, '%f %f %f\n', X{n}');
end
fclose(fid);
