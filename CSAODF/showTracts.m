% Visualizes the computed tracts. An optional parameter will be ignored if
% [] is entered for it.
%
% showTracts(tracts, overlayVol, nTracts, mask, threshScore, showSeeds)
%  
% tracts:       Tract structure resulting from HoughTract.m.
% overlayVol:   Optional 3D image over which the tracts will be visualized.
% nTracts:      Optional number of top-scored tracts to be shown.
% mask:         Optional mask(s), either a single 3D image or a cell of
%               multiple 3D images, each containing a region of interest.
%               If provided, only the tracts that intersect all of the
%               masks will be shown.
% threshScore:  If provided, only the tracts with a score above
%               'threshScore' will be shown.
% showSeeds:    If true, seed points will be drawn as circles.
%
% See also:  HoughTract, tracts2vox, exportTrackVis, exportText, showODFs,
%            connMatrix, CSAODF_CLI, EXAMPLE, EXAMPLE_CLI.

% Code by Iman Aganj.

function p = showTracts(tracts, overlayVol, nTracts, mask, threshScore, showSeeds)

score = tracts.score;
if exist('mask', 'var') && ~isempty(mask)
    [X, seedPointsVox, ind, mask] = tracts2vox(tracts, mask);
    score = score(ind);
else
    [X, seedPointsVox] = tracts2vox(tracts);
end
if ~exist('nTracts', 'var') || isempty(nTracts)
    nTracts = length(X);
else
    nTracts = min(nTracts, length(X));
end
[score, ind] = sort(score, 'descend');
ind = ind(1:nTracts);
score = score(1:nTracts);
X = X(ind);
seedPointsVox = seedPointsVox(ind, :);
if ~exist('threshScore', 'var') || isempty(threshScore)
    threshScore = min(score(~isinf(score)));
end
nJet = nTracts*100;
map = jet(nJet);
score = score - max(min(score(:)), threshScore); score = score/max(score(:));

figure
hold on
axis equal tight
set(gca, 'Clipping', 'off')
xlabel X, ylabel Y, zlabel Z
if exist('overlayVol', 'var') && ~isempty(overlayVol)
    overlayVol = permute(overlayVol, [2 1 3]);
    centerVol = round((size(overlayVol)+1)/2);
    slice(overlayVol, centerVol(1), centerVol(2), centerVol(3));
    shading interp
    colormap gray
    alpha 0.3
end
nShownTracts = 0;
for i = 1:nTracts
    if score(i)>=0 && size(X{i},1)>1
        cl = round((nJet-1)*score(i)) + 1;
        if exist('showSeeds', 'var') && showSeeds
            plot3(seedPointsVox(i,1), seedPointsVox(i,2), seedPointsVox(i,3), 'o', 'color', map(cl,:))
        end
        if nargout > 0
            p(ind(i),1) = plot3(X{i}(:,1), X{i}(:,2), X{i}(:,3), 'color', map(cl,:));
        else
            plot3(X{i}(:,1), X{i}(:,2), X{i}(:,3), 'color', map(cl,:));
        end
        nShownTracts = nShownTracts + 1;
    end
    if mod(i,1000)==0
        drawnow
    end
end
if exist('mask', 'var') && ~isempty(mask)
    nMasks = length(mask);
    for i=1:nMasks
        if nMasks==1
            ratio = .5;
        else
            ratio = (i-1)/(nMasks-1);
        end
        patch(isosurface(permute(mask{i}, [2 1 3]), .5), 'FaceAlpha', .4, 'EdgeColor', 'none', 'FaceColor', map(round((nJet-1)*ratio)+1,:));
    end
    camlight
    lighting gouraud
end
hold off
set(gcf, 'Name', [num2str(nShownTracts) ' tracts'])
