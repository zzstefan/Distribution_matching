% Augments the connectivity matrix with indirect connections by modeling
% connectivity as electric conductance.
%
% aC = augConnMatrix(C)
%
% C:         Connectivity matrix (symmetric).
%
% aC:        Augmented connectivity matrix.
%
% See also:  connMatrix, FreeSurferLUT, HoughTract, EXAMPLE.

% Code by Iman Aganj.
%
% Reference:
% I. Aganj, G. Prasad, P. Srinivasan, A. Yendiki, P. M. Thompson, and B.
% Fischl, "Structural brain network augmentation via Kirchhoff's laws,"
% in Proc. Joint Annual Meeting of ISMRM-ESMRMB, Milan, Italy, 2014.
% https://cds.ismrm.org/protected/14MProceedings/PDFfiles/2665.pdf
%
% A more direct method for conductance-based quantification of brain
% connectivity is available at: www.nitrc.org/projects/conductance

function aC = augConnMatrix(C)

N = size(C,1);
if size(C,2)~=N
    error('C has to be a square matrix!')
end

% Laplacian matrix:
C(logical(eye(N))) = 0;
L = diag(sum(C,2)) - C;
L(1,:) = [1 zeros(1,N-1)];  % Node 1 is grounded.
invL = inv(L);

% Computing the total conductance matrix:
aC = inf(N);
for i = 1:N
    for j = setdiff(1:N,i)
        I = zeros(N,1);
        I([i j]) = [1, -1];
        aC(i,j) = 1/(I'*invL*[0; I(2:end)]);  % Node 1 is grounded.
    end
end

if any(isnan(aC(:)))
    warning('NaN values were produced, possibly because the input matrix has disconnected components. Adding 1e-10 to ''C'' and recomputing...')
    aC = augConnMatrix(C + 1e-10);
end

aC(logical(eye(N))) = nan;
