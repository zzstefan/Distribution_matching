% Computes the q-ball CSA-ODFs from a single shell of diffusion MRI data.
%
% [sh, mask] = reconCSAODF(S, S0, angles, param, mask, pT)
%  
% S:         4D array (or 1D cell array of 3D images) containing the
%            diffusion-weighted images. 
% S0:        Non-diffusion-weighted (b=0) image. Use [] if S0=1 everywhere.
% angles:    A matrix with two columns [theta,phi] of the diffusion
%            directions. (cart2sph_phys.m can be used to convert the
%            gradient table [x,y,z] to [theta,phi].)
% param:     Optional structure containing the following parameters.
%            A default value is set to any non-existing field.
%    param.method:       Reconstruction method
%                           1 --> CSA-ODF (default).
%                           2 --> Original (Tuch's) ODF.
%                           3 --> Original (Tuch's) ODF, with
%                                 Laplace-Beltrami sharpening (lambda=0.1).
%    param.basisOrder:   Spherical harmonic basis order (an even number,
%                        default: 4). 
%    param.delta:        Regularization parameter (default: 0.0001).
%    param.indexGPU:     If provided, the GPU device specified by this
%                        index will be used. Run 'CSAODF_CLI GPUs' to see
%                        the indices of the available GPUs.
% mask:      Optional image mask (can be ignore by entering []). It will be
%            modified in the output to contain only voxels with S0>0. 
% pT:        This matrix can be pre-calculated and optionally provided to
%            speed up the computation (see the code for details). 
% 
% sh:        Spherical harmonic coefficients of the ODFs.
%
% See also:  cart2sph_phys, reconCSAODF3Q, reconDT, showODFs, sampleODFs,
%            makeGFA, HoughTract, CSAODF_CLI, EXAMPLE, EXAMPLE_CLI.

% Codes by Iman Aganj.
%
% Reference:
% I. Aganj, C. Lenglet, G. Sapiro, E. Yacoub, K. Ugurbil, and N. Harel,
% "Reconstruction of the orientation distribution function in single and
% multiple shell q-ball imaging within constant solid angle,"
% Magnetic Resonance in Medicine, vol. 64, no. 2, pp. 554-566, 2010.

function [sh, mask] = reconCSAODF(S, S0, angles, param, mask, pT)

if ~exist('param', 'var')
    param = [];
end
if ~isfield(param, 'basisOrder')
    param.basisOrder = 4;
end
if ~isfield(param, 'delta')
    param.delta = 1e-4;
end
if ~isfield(param, 'method')
    param.method = 1;
end
if ~isfield(param, 'indexGPU')
    param.indexGPU = 0;
end
existS0 = ~isempty(S0);
if existS0
    S0 = double(S0);
end
if param.indexGPU
    gpu = gpuDevice(param.indexGPU);
    disp(['Using the "' gpu.Name '" GPU...']);
    if existS0
        S0 = gpuArray(S0);
    end
    delta = gpuArray(param.delta);
else
    delta = param.delta;
end
if iscell(S)
    sz = size(S{1});
    if length(sz)==2
        sz = [sz 1];
    end
    nGrads = length(S);
    if param.indexGPU
        S = cellfun(@gpuArray, S, 'UniformOutput', false);
    end
else
    S = double(S);
    if param.indexGPU
        S = gpuArray(S);
    end
    if ndims(S)<4
        warning('S must be a 4D array. Shifting the 3rd dimension...')
        S = permute(S, [1 2 4 3]);
    end
    sz = size(S);
    nGrads = sz(4);
    sz = sz(1:3);
end
if size(angles,1) ~= nGrads
    error('''angles'' should have the same number of rows as the number of diffusion images in ''S''!')
end
if existS0
    if exist('mask', 'var') && ~isempty(mask)
        mask = mask & S0>0;
    else
        mask = S0>0;
    end
end
[L, C] = makeL(param.basisOrder);
nSH = size(L,1);
if ~exist('pT', 'var') || isempty(pT) % pT can be precomputed and provided as input to gain speed.
    [~, pT] = makePT(param.basisOrder, angles);
    %[~, pT] = makePT(param.basisOrder, angles, .0000005, L); % Use this instead to incorporate the regularization in (Descoteaux et al, MRM'07).
end
switch param.method
    case 1
        iT = diag(C.*L/(16*pi^2)) * pT;
    case 2
        iT = diag(C) * pT;
    case 3
        iT = diag((1-.1*L).*C) * pT;
    otherwise
        error('Choose 1, 2, or 3 for method.')
end

if iscell(S)
    sampleS = S{1}(1);
else
    sampleS = S(1);
end
sh = zeros([sz nSH], 'like', sampleS);
matRow = zeros(1,1,1,nSH, 'like', sampleS);
oldVer = verLessThan('matlab', '9.1'); % Matlab releases earlier than R2016b.
for n = 1:nGrads
    if existS0
        if iscell(S)
            rE = reg(double(S{n})./S0, delta);
        else
            rE = reg(S(:,:,:,n)./S0,   delta);
        end
    else
        if iscell(S)
            rE = reg(double(S{n}), delta);
        else
            rE = reg(S(:,:,:,n),   delta);
        end
    end
    if param.method == 1
        rE = log(-log(rE));
    end
    matRow(:) = iT(:,n);
    if oldVer
        sh = sh + bsxfun(@times, rE, matRow);
    else
        sh = sh + rE .* matRow;
    end
end
if exist('mask', 'var') && ~isempty(mask)
    if oldVer
        sh(:,:,:,2:end) = bsxfun(@times, sh(:,:,:,2:end), mask);
    else
        sh(:,:,:,2:end) = sh(:,:,:,2:end) .* mask;
    end
else
    mask = [];
end
if param.method > 1
    sh = bsxfun(@rdivide, sh, sh(:,:,:,1)) / (2*sqrt(pi));
end
sh(isnan(sh)) = 0;
sh(:,:,:,1) = 1/(2*sqrt(pi));
if param.indexGPU
    sh = gather(sh);
    mask = gather(mask);
end
