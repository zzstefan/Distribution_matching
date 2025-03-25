% Creates the diagonal entries of the Laplacian 'L' and the Funk-Radon
% Transform 'C' matrices for the real and symmetric spherical harmonic
% basis of order 'basisOrder'.
% 
% [L, C] = makeL(basisOrder)
%
% See also:  reconCSAODF, reconCSAODF3Q, CSAODF_CLI, EXAMPLE, EXAMPLE_CLI.
 
% Code by Iman Aganj.

function [L, C] = makeL(basisOrder)

L = zeros((basisOrder+1)*(basisOrder+2)/2, 1);
C = L;
for k = 0:2:basisOrder
    for m = -k:k
        j = k*(k+1)/2 + m + 1;
        L(j) = -k*(k+1);
        C(j) = ((-1)^(k/2))*prod(1:2:(k-1))/prod(2:2:k);
    end
end
C = (2*pi)*C;
