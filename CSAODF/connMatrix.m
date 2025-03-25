% Creates the connectivity matrix from tractography results and labels, by
% counting the number (weighted by score) of tracts passing through every
% pair of ROIs. Includes interactive visualization.
%
% [C, labels] = connMatrix(tracts, labelVol, labels, weightByScore, showConnTracts, labelNames)
%
% tracts:         Tract structure resulting from HoughTract.m.
% labelVol:       Segmentation volume, including all the labels.
% labels (in):    Optional vector of label numbers to be included in the
%                 connectivity matrix (default: all labels except 0).
% weightByScore:  If true (default), tracts will be weighted by their scores.
% showConnTracts: If true, the connectivity matrix and the tracts will be
%                 visualized (default: false). By moving the mouse pointer
%                 on the matrix, the tractogram will be updated to show
%                 only the tracts connecting the corresponding pair of
%                 ROIs. Clicking will lock the pointer's position.
% labelNames:     An optional cell array with the same number of elements
%                 as in 'labels', where each cell contains the name of the
%                 corresponding label. If 'showConnTracts' is true and
%                 'labelNames' is provided, the label names will be used in
%                 the visualization.
%
% C:              Symmetric connectivity matrix.
% labels (out):   Labels included in the matrix (same as input).
%
% See also:  augConnMatrix, FreeSurferLUT, HoughTract, showTracts,
%            exportTrackVis, exportText, tracts2vox, EXAMPLE.
 
% Code by Iman Aganj.

function [C, labels] = connMatrix(tracts, labelVol, labels, weightByScore, showConnTracts, labelNames)

if ~exist('labels', 'var') || isempty(labels)
    labels = unique(labelVol(:));
    if ismember(0, labels)
        labels = setdiff(labels, 0);
        warning('Label 0 was excluded.')
    end
end
if ~exist('weightByScore', 'var') || isempty(weightByScore)
    weightByScore = true;
end
nTracts = length(tracts.score);
nROIs = length(labels);

X = tracts2vox(tracts);
indM1 = cellfun(@(x) size(x,1)>1, X);
X(indM1) = cellfun(@(x) interp1(x,1:.2:size(x,1)), X(indM1), 'UniformOutput', false); % Upsample tracts 5x

interROI = zeros(nTracts, nROIs);
for i = 1:nROIs
    [~, ~, ind] = tracts2vox(X, {labelVol==labels(i)});
    interROI(ind,i) = 1;
end
clear X indM1
if weightByScore
    C = bsxfun(@times, interROI, subplus(tracts.score))' * interROI;
else
    C = interROI' * interROI;
end
C(logical(eye(nROIs))) = nan;

if exist('showConnTracts', 'var') && showConnTracts
    [~,~,ind] = unique(labelVol);
    overlayVol = zeros(size(labelVol));
    overlayVol(:) = ind;
    hLine = showTracts(tracts, overlayVol);
    isLine = arrayfun(@(l) isa(l, 'matlab.graphics.chart.primitive.Line'), hLine);
    hold on
    hPatch = arrayfun(@(l) patch(isosurface(permute(labelVol==l, [2 1 3]), .5), 'FaceAlpha', .4, 'EdgeColor', 'none', 'Visible', false), labels);
    hold off
    set(gcf, 'units', 'normalized')
    set(gcf, 'Position', [0 0 .5 1])
    
    figure('units', 'normalized', 'Position', [.5 0 .5 1], 'WindowButtonMotionFcn', @mouseMove, 'WindowButtonDownFcn', @mouseButton, 'Name', 'Move your mouse to mask the tracts. Click to lock.');
    imagesc(C)
    axis equal tight
    colormap(jet(2^12))
    colorbar
    ax = gca;
    ax.UserData.locked = false;
    if exist('labelNames', 'var')
        set(ax, 'TickLabelInterpreter', 'none')
        xticks(1:nROIs)
        xticklabels(labelNames)
        yticks(1:nROIs)
        yticklabels(labelNames)
    end
end

    function mouseMove(hObject, eventData)
        if ~ax.UserData.locked
            x = round(hObject.CurrentAxes.CurrentPoint(1, [2 1]));
            x = min(max(x,1),nROIs);
            interLines = all(interROI(:,x),2);
            set(hLine(isLine & ~interLines), 'Visible', false)
            set(hLine(isLine & interLines), 'Visible', true)
            set(hPatch(x(1)), 'FaceColor', 'r')
            set(hPatch(x(2)), 'FaceColor', 'b')
            set(hPatch(setdiff(1:nROIs,x)), 'Visible', false)
            set(hPatch(x), 'Visible', true)
            if exist('labelNames', 'var')
                ax.Title.String = ['\color{red}' labelNames{x(1)} ' \color{black}- \color{blue}' labelNames{x(2)}];
            else
                ax.Title.String = ['\color{red}' num2str(x(1)) ' \color{black}- \color{blue}' num2str(x(2))];
            end
        end
    end

    function mouseButton(hObject, eventData)
        ax.UserData.locked = ~ax.UserData.locked;
        if ax.UserData.locked
            ax.XLabel.String = 'Locked';
        else
            ax.XLabel.String = '';
            mouseMove(hObject)
        end
    end

end
