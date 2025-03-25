% Computes the diffusion tensors from diffusion MRI data.
%
% [D, FA, MD, mask, sh] = reconDT(S, S0, angles, b, param, mask)
%  
% S:         4D array (or 1D cell array of 3D images) containing the
%            diffusion-weighted images. 
% S0:        Non-diffusion-weighted (b=0) image. Use [] if S0=1 everywhere.
% angles:    A matrix with two columns [theta,phi] of the diffusion
%            directions. (cart2sph_phys.m can be used to convert the
%            gradient table [x,y,z] to [theta,phi].)
% b:         B-value, either a scalar or a vector with values corresponding
%            to rows of angles.
% param:     Optional structure containing the following parameters.
%            A default value is set to any non-existing field.
%    param.delta:        Regularization parameter (default: 0.0001).
%    param.indexGPU:     If provided, the GPU device specified by this
%                        index will be used. Run 'CSAODF_CLI GPUs' to see
%                        the indices of the available GPUs.
% mask:      Optional image mask (can be ignore by entering []). It will be
%            modified in the output to contain only voxels with S0>0. 
% 
% D:         Diffusion tensor coefficients.
% FA:        Fractional anisotropy.
% MD:        Mean diffusivity.
% sh:        Spherical harmonic coefficients of ODFs made from the tensors.
%
% See also:  cart2sph_phys, DT2ODF, reconCSAODF, reconCSAODF3Q.

% Codes by Iman Aganj.

function [D, FA, MD, mask, sh] = reconDT(S, S0, angles, b, param, mask)

if ~exist('param', 'var')
    param = [];
end
if ~isfield(param, 'delta')
    param.delta = 1e-4;
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
    b = gpuArray(b);
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
if isscalar(b)
    b = b * ones(nGrads,1);
elseif length(b) ~= nGrads
    error('''b'' should be either a scalar or a vector with the same number of elements as the number of diffusion images in ''S''!')
end
if existS0
    if exist('mask', 'var') && ~isempty(mask)
        mask = mask & S0>0;
    else
        mask = S0>0;
    end
end

bvecs = cart2sph_phys(angles, true);
T = [bvecs(:,1).^2  2*bvecs(:,1).*bvecs(:,2)  2*bvecs(:,1).*bvecs(:,3)  bvecs(:,2).^2  2*bvecs(:,2).*bvecs(:,3)  bvecs(:,3).^2];
iT = pinv(T);

if iscell(S)
    sampleS = S{1}(1);
else
    sampleS = S(1);
end
D = zeros([sz 6], 'like', sampleS);
matRow = zeros(1,1,1,6, 'like', sampleS);
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
    rE = -log(rE) / b(n);
    matRow(:) = iT(:,n);
    if oldVer
        D = D + bsxfun(@times, rE, matRow);
    else
        D = D + rE .* matRow;
    end
end
clear S S0
if exist('mask', 'var') && ~isempty(mask)
    if oldVer
        D = bsxfun(@times, D, mask);
    else
        D = D .* mask;
    end
else
    mask = [];
end
D(isnan(D)) = 0;
if param.indexGPU
    D = gather(D);
    mask = gather(mask);
end
if nargout > 1
    FA = real(sqrt(1.5 - .5*(sum(D(:,:,:,[1 4 6]),4).^2)./sum(D.^2,4)));
    if nargout > 2
        MD = mean(D(:,:,:,[1 4 6]),4);
        if nargout > 4
            sh = DT2ODF(D, param);
        end
    end
end
