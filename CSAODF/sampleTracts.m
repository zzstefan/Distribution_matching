% Makes a new tract structure from 'tracts' containing only the tracts with
% indices specified in 'indSeeds' (which can be generated using tracts2vox.m).
%
% tracts = sampleTracts(tracts, indSeeds)
%  
% See also:  HoughTract, tracts2vox, showTracts, exportTrackVis,
%            exportText, connMatrix, CSAODF_CLI, EXAMPLE, EXAMPLE_CLI.

% Code by Iman Aganj.

function tracts = sampleTracts(tracts, indSeeds)

tracts.seedPoints = tracts.seedPoints(indSeeds, :);
tracts.partialLengths = tracts.partialLengths(indSeeds, :);
tracts.polyCoeffs= tracts.polyCoeffs(indSeeds, :);
tracts.score = tracts.score(indSeeds);
