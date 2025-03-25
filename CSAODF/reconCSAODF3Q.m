% Computes the q-ball CSA-ODFs from three shells of diffusion MRI data. The
% b-values of the shells must be such that 0,b1,b2,b3 is an arithmetic
% sequence, and the gradient orientations must be the same for all three
% shells.
%
% [sh, mask] = reconCSAODF3Q(S1, S2, S3, S0, angles, param, mask, pT)
%  
% S1,S2,S3:  4D arrays (or 1D cell arrays of 3D images) containing the
%            diffusion-weighted images of each shell. 
% S0:        Non-diffusion-weighted (b=0) image.
% angles:    A matrix with two columns [theta,phi] of the diffusion
%            directions. (cart2sph_phys.m can be used to convert the
%            gradient table [x,y,z] to [theta,phi].)
% param:     Optional structure containing the following parameters. A
%            default value is set to any non-existing field. 
%    param.basisOrder:   Spherical harmonic basis order (an even number, default: 4).
%    param.delta:        Regularization parameter (default: 0.0001).
%    param.delta3:       Projection regularization parameter (default: 0.01).
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
% See also:  cart2sph_phys, reconCSAODF, reconDT, showODFs, sampleODFs,
%            makeGFA, HoughTract, CSAODF_CLI, EXAMPLE, EXAMPLE_CLI.

% Codes by Iman Aganj.
%
% Reference:
% I. Aganj, C. Lenglet, G. Sapiro, E. Yacoub, K. Ugurbil, and N. Harel,
% "Reconstruction of the orientation distribution function in single and
% multiple shell q-ball imaging within constant solid angle,"
% Magnetic Resonance in Medicine, vol. 64, no. 2, pp. 554-566, 2010.

function [sh, mask] = reconCSAODF3Q(S1, S2, S3, S0, angles, param, mask, pT)

if ~exist('param', 'var')
    param = [];
end
if ~isfield(param, 'basisOrder')
    param.basisOrder = 4;
end
if ~isfield(param, 'delta')
    param.delta = 1e-4;
end
if ~isfield(param, 'delta3')
    param.delta3 = 0.01;
end
if ~isfield(param, 'indexGPU')
    param.indexGPU = 0;
end
sz = size(S0);
if length(sz)==2
    sz = [sz 1];
end
S0 = double(S0);
if param.indexGPU
    gpu = gpuDevice(param.indexGPU);
    disp(['Using the "' gpu.Name '" GPU...']);
    S0 = gpuArray(S0);
    delta = gpuArray(param.delta);
    delta3 = gpuArray(param.delta3);
    sqrt5 = gpuArray(sqrt(5));
    sqrt6 = gpuArray(sqrt(6));
else
    delta = param.delta;
    delta3 = param.delta3;
    sqrt5 = sqrt(5);
    sqrt6 = sqrt(6);
end
if iscell(S1)
    nGrads = length(S1);
    if param.indexGPU
        S1 = cellfun(@gpuArray, S1, 'UniformOutput', false);
        S2 = cellfun(@gpuArray, S2, 'UniformOutput', false);
        S3 = cellfun(@gpuArray, S3, 'UniformOutput', false);
    end
else
    S1 = double(S1);
    S2 = double(S2);
    S3 = double(S3);
    if param.indexGPU
        S1 = gpuArray(S1);
        S2 = gpuArray(S2);
        S3 = gpuArray(S3);
    end
    if ndims(S1)<4
        warning('S1, S2, and S3 must be 4D arrays. Shifting the 3rd dimension...')
        S1 = permute(S1, [1 2 4 3]);
        S2 = permute(S2, [1 2 4 3]);
        S3 = permute(S3, [1 2 4 3]);
    end
    nGrads = size(S1, 4);
end
if size(angles,1) ~= nGrads
    error('''angles'' should have the same number of rows as the number of diffusion images in ''S''!')
end
if exist('mask', 'var') && ~isempty(mask)
    mask = mask & S0>0;
else
    mask = S0>0;
end
oldVer = verLessThan('matlab', '9.1'); % Matlab releases earlier than R2016b.
[L, C] = makeL(param.basisOrder);
nSH = size(L,1);
if ~exist('pT', 'var') || isempty(pT) % pT can be precomputed and provided as input to gain speed.
    [~, pT] = makePT(param.basisOrder, angles);
    %[~, pT] = makePT(param.basisOrder, angles, .0000005, L); % Use this instead to incorporate the regularization in (Descoteaux et al, MRM'07).
end
iT = diag(C.*L/(16*pi^2)) * pT;
sh = zeros([sz nSH], 'like', S0);
matRow = zeros(1,1,1,nSH, 'like', S0);
deltaVol = delta * ones(sz);
for n = 1:nGrads
    if iscell(S1)
        E1 = -log(reg(double(S1{n})./S0, delta));
        E2 = -log(reg(double(S2{n})./S0, delta));
        E3 = -log(reg(double(S3{n})./S0, delta));
    else
        E1 = -log(reg(S1(:,:,:,n)./S0, delta));
        E2 = -log(reg(S2(:,:,:,n)./S0, delta));
        E3 = -log(reg(S3(:,:,:,n)./S0, delta));
    end
    % Projecting the values of E1,E2,E3 onto the subspace defined by inequalities [31], and computing the parameters according to Eqs. [29-30].
    aux1 = E1 + E2 + E3;
    aux2 = E2 - 2*E1;
    aux3 = 2*E2 + E1;
    C01 = (1/6+0.5*(sqrt5+1)*delta3)*aux1<E1;
    C02 = E1<(1/3-(sqrt5+2)*delta3)*aux1;
    C03 = aux2<=(-1/3+(2*sqrt5+5)*delta3)*aux1;
    C04 = aux2<(-sqrt5*delta3)*aux1;
    C05 = aux3<(1-sqrt5*delta3)*aux1;
    C06 = (5/6+0.5*(5+sqrt5)*delta3)*aux1<aux3;
    C1 = E2<(1/3+delta3)*aux1 & C01 & C02;
    C2 = C03 & ~C02;
    C3 = aux3>(1-sqrt5*delta3)*aux1 & ~C03 & C04;
    C4 = ~C05 & ~C04;
    C5 = C06 & C05 & aux2>(-sqrt5*delta3)*aux1;
    C6 = ~C01 & ~C06;
    C7 = ~(C1 | C2 | C3 | C4 | C5 | C6);
    E1 = exp(-(E1.*C1 + (1/3-(sqrt5+2)*delta3)*aux1.*C2 + ((.2-delta3/sqrt5)*aux1-.4*aux2).*C3 + (.2+delta3/sqrt5)*aux1.*C4 + (.2*aux3+(2*delta3/sqrt5)*aux1).*C5 + (1/6+.5*(sqrt5+1)*delta3)*aux1.*C6 + E1.*C7));
    E2 = exp(-((1/3+delta3)*aux1.*C1 + (1/3+delta3)*aux1.*C2 + ((.4-2*delta3/sqrt5)*aux1+.8*aux2).*C3 + (.4-3*delta3/sqrt5)*aux1.*C4 + (.4*aux3-(delta3/sqrt5)*aux1).*C5 + (1/3+delta3)*aux1.*C6 + E2.*C7));
    E3 = exp(-aux1)./(E1.*E2);
    aux1 = E2-E1.^2;
    aux2 = (E3-E1.*E2)./(2*aux1);
    aux3 = aux2.^2 - (E1.*E3-E2.^2)./aux1;
    aux1 = sqrt(aux1.*(aux1>0));
    aux3 = sqrt(aux3.*(aux3>0));
    alpha = aux2 + aux3;
    beta  = aux2 - aux3;
    aux14 = cat(4, aux1, aux1, aux1, deltaVol, (aux1+alpha-beta-sqrt6*deltaVol)/3, deltaVol, deltaVol, deltaVol, aux1, .2*(2*alpha+aux1-2*(sqrt6+1)*deltaVol), .2*(-2*beta+aux1+2-2*(sqrt6+1)*deltaVol), deltaVol, deltaVol, deltaVol, .5-(1+sqrt6/2)*deltaVol);
    alpha4 = cat(4, alpha, alpha, 1-deltaVol, alpha, (2*aux1+5*alpha+beta+sqrt6*deltaVol)/6, alpha, 1-deltaVol, .5*(alpha+beta)+(1+sqrt6/2)*deltaVol, 1-deltaVol, .2*(4*alpha+2*aux1+(sqrt6+1)*deltaVol), 1-deltaVol, (sqrt6+3)*deltaVol, 1-deltaVol, 1-deltaVol, 1-deltaVol);
    beta4 = cat(4, beta, deltaVol, beta, beta, (-2*aux1+alpha+5*beta-sqrt6*deltaVol)/6, deltaVol, beta, .5*(alpha+beta)-(1+sqrt6/2)*deltaVol, deltaVol, deltaVol, .2*(4*beta-2*aux1+1-(sqrt6+1)*deltaVol), deltaVol, deltaVol, 1-(sqrt6+3)*deltaVol, deltaVol);
    if oldVer
        dist2 = bsxfun(@minus, aux14, aux1).^2 + bsxfun(@minus, alpha4, alpha).^2 + bsxfun(@minus, beta4, beta).^2;
    else
        dist2 = (aux14 - aux1).^2 + (alpha4 - alpha).^2 + (beta4 - beta).^2;
    end
    dist2(~(aux14>.99*delta & alpha4-beta4-2*aux14>.99*delta*sqrt6 & beta4>.99*delta & alpha4<1-.99*delta)) = inf;
    [~, ind] = min(dist2, [], 4);
    ind = (1:prod(sz))' + (ind(:)-1)*prod(sz);
    aux1(:) = aux14(ind);
    alpha(:) = alpha4(ind);
    beta(:) = beta4(ind);
    lambda = .5*(1 + sqrt(1-(2*aux1./(alpha-beta)).^2));
    C01 = abs(lambda.*(alpha-beta)+beta-E1) + abs(lambda.*(alpha.^2-beta.^2)+beta.^2-E2) + abs(lambda.*(alpha.^3-beta.^3)+beta.^3-E3) < abs((1-lambda).*(alpha-beta)+beta-E1) + abs((1-lambda).*(alpha.^2-beta.^2)+beta.^2-E2) + abs((1-lambda).*(alpha.^3-beta.^3)+beta.^3-E3);
    lambda = lambda.*C01 + (1-lambda).*~C01;
    matRow(:) = iT(:,n);
    if oldVer
        sh = sh + bsxfun(@times, lambda.*log(-log(alpha)) + (1-lambda).*log(-log(beta)), matRow);
    else
        sh = sh + (lambda.*log(-log(alpha)) + (1-lambda).*log(-log(beta))) .* matRow;
    end
end
sh(:,:,:,1) = 1/(2*sqrt(pi));
if oldVer
    sh(:,:,:,2:end) = bsxfun(@times, sh(:,:,:,2:end), mask);
else
    sh(:,:,:,2:end) = sh(:,:,:,2:end) .* mask;
end
sh(isnan(sh)) = 0;
if param.indexGPU
    sh = gather(sh);
    mask = gather(mask);
end
