% Computes the generalized fractional anisotropy 'GFA' from the 4D array of
% ODFs 'sh'.
%
% GFA = makeGFA(sh)
%  
% See also:  reconCSAODF, reconCSAODF3Q, HoughTract, CSAODF_CLI, EXAMPLE,
%            EXAMPLE_CLI.

% Code by Iman Aganj.

function GFA = makeGFA(sh)

GFA = 1./sqrt(1+1./(4*pi*sum(sh(:,:,:,2:end).^2,4)));
