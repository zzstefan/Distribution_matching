% Generates a volume from the computed 'tracts', with each voxel indicating
% the number of tracts that pass through it. If the optional
% 'weightedByScore' is true, the sum will be weighted by the tract scores.
%
% V = fiberDensityVol(tracts, weightedByScore)
%  
% See also:  HoughTract, tracts2vox, exportTrackVis, exportText,
%            showTracts, connMatrix, CSAODF_CLI, EXAMPLE.

% Code by Iman Aganj.

function V = fiberDensityVol(tracts, weightedByScore)

if ~exist('weightedByScore', 'var')
    weightedByScore = false;
end
X = tracts2vox(tracts);
sz = tracts.param.sz(1:3);
V = zeros(sz);
for n = 1:length(X)
    x = int32(X{n});
    legitX = all(x>0, 2) & all(bsxfun(@le, x, sz),2);
    ind = sub2ind(sz, x(legitX,1), x(legitX,2), x(legitX,3));
    if weightedByScore
        V(ind) = V(ind) + tracts.score(n);
    else
        V(ind) = V(ind) + 1;
    end
end
