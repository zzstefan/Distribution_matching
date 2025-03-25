% Samples the 4D ODF array 'sh' along the directions specified by two
% columns (theta,phi) in 'angles', and returns the 4D array 'odf' of
% the sampled values.
% 
% odf = sampleODFs(sh, angles)
%
% See also:  compSH, cart2sph_phys, reconCSAODF, reconCSAODF3Q, showODFs,
%            CSAODF_CLI, EXAMPLE, EXAMPLE_CLI.

% Code by Iman Aganj.

function odf = sampleODFs(sh, angles)

if verLessThan('matlab', '9.1') % Matlab releases earlier than R2016b.
    odf = sum(bsxfun(@times, permute(sh, [1 2 3 5 4]), permute(compSH((sqrt(8*size(sh,4)+1)-3)/2, angles(:,1)', angles(:,2)'), [3 4 5 2 1])), 5);
else
    odf = sum(permute(sh, [1 2 3 5 4]) .* permute(compSH((sqrt(8*size(sh,4)+1)-3)/2, angles(:,1)', angles(:,2)'), [3 4 5 2 1]), 5);
end

