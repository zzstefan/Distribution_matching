% Uses the FreeSurfer lookup table to find label name(s).
%
% labelName = FreeSurferLUT(N)
%
% N:          An array of FreeSurfer label numbers.
%
% labelName:  A cell array of label names (or just the label name if 'N' is
%             a scalar).
%
% See also:   connMatrix, augConnMatrix, tracts2vox, showTracts,
%             exportTrackVis, exportText, EXAMPLE.
 
% Code by Iman Aganj.

function labelName = FreeSurferLUT(N, Cbm)

if exist('Cbm', 'var') && Cbm
    fid = fopen('CerebellumColorLUT.txt');
else
    fid = fopen('FreeSurferColorLUT.txt');
end
S = textscan(fid, '%d %s %*[^\n]', 'HeaderLines', 3, 'CommentStyle', '#');
labelName = arrayfun(@(n) S{2}{find(S{1}==n,1)}, N, 'UniformOutput', false);

if length(labelName)==1
    labelName = labelName{1};
end

fclose(fid);
