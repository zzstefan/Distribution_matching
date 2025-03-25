% Computes the CSA-ODFs from diffusion tensors.
%
% sh = DT2ODF(D, param, res)
%  
% D:         4D array of diffusion tensor coefficients.
% param:     Optional structure containing the following parameters.
%            A default value is set to any non-existing field.
%    param.method:       Reconstruction method
%                           1 --> CSA-ODF (default).
%                           2 --> Original (Tuch's) ODF.
%                           3 --> Original (Tuch's) ODF, with
%                                 Laplace-Beltrami sharpening (lambda=0.1).
%    param.basisOrder:   Spherical harmonic basis order (an even number,
%                        default: 2). 
%    param.delta:        Regularization parameter (default: 0.0001).
%    param.indexGPU:     If provided, the GPU device specified by this
%                        index will be used. Run 'CSAODF_CLI GPUs' to see
%                        the indices of the available GPUs.
% res:       Optional angular resolution in radians (default: 0.5).
% 
% sh:        Spherical harmonic coefficients of the ODFs.
%
% See also:  reconDT, cart2sph_phys, reconCSAODF, reconCSAODF3Q, showODFs,
%            sampleODFs, makeGFA, HoughTract, CSAODF_CLI, EXAMPLE, EXAMPLE_CLI.
 
% Codes by Iman Aganj.
%
% Reference:
% I. Aganj, C. Lenglet, G. Sapiro, E. Yacoub, K. Ugurbil, and N. Harel,
% "Reconstruction of the orientation distribution function in single and
% multiple shell q-ball imaging within constant solid angle,"
% Magnetic Resonance in Medicine, vol. 64, no. 2, pp. 554-566, 2010.

function sh = DT2ODF(D, param, res)

if ~exist('param', 'var')
    param = [];
end
if ~exist('res', 'var') || isempty(res)
    res = .5;
end
if ~isfield(param, 'basisOrder')
    param.basisOrder = 2;
end
if ~isfield(param, 'indexGPU')
    param.indexGPU = 0;
end
if ndims(D)<4
    warning('D must be a 4D array. Shifting the 3rd dimension...')
    D = permute(D, [1 2 4 3]);
end
if size(D,4) ~= 6
    error('The fourth dimension of ''D'' must be 6!')
end
b = 1000;
nDir = 0;
for theta = res:res:(pi-res)
    for phi = 0:(res/sin(theta)):2*pi
        nDir = nDir + 1;
        sampleAngles(nDir,:) = [theta phi];
    end
end
bvecs = cart2sph_phys(sampleAngles, true);
if param.indexGPU
    gpuDevice(param.indexGPU);
    D = gpuArray(D);
    bvecs = gpuArray(bvecs);
    b = gpuArray(b);
end

T = [bvecs(:,1).^2  2*bvecs(:,1).*bvecs(:,2)  2*bvecs(:,1).*bvecs(:,3)  bvecs(:,2).^2  2*bvecs(:,2).*bvecs(:,3)  bvecs(:,3).^2];
T = permute(T, [3 4 5 1 2]);

S = 0;
for i = 1:6
    S = S + D(:,:,:,i) .* T(1,1,1,:,i);
end
clear D
S = exp(-b*S);
sh = reconCSAODF(S, [], sampleAngles, param);
