% Verifies if the gradient table uses the same coordinates as the image.
% An optional parameter will be ignored if [] is entered for it.
%
% [optPerm, optFlip, grad, allErrors] = checkGradTable(S, S0, grad, Mask, indexGPU, voxelSize, res, reconParam)
%  
% S:          4D array containing the diffusion-weighted images. 
% S0:         Non-diffusion-weighted (b=0) image. Use [] if S0=1 everywhere.
% grad (in):  Gradient table (tall matrix).
% Mask:       Optional fibrous tissue mask. If not specified (not
%             recommended), a rough mask will be created, in which case,
%             the b-value must be specified in reconParam.bValue.
% indexGPU:   If provided, the GPU device specified by this index will be
%             used. Run 'CSAODF_CLI GPUs' to see the indices of the
%             available GPUs.
% voxelSize:  Optional vector of the voxel size (default: [1 1 1]).
% res:        Optional angular resolution in radians (default: 0.5).
% reconParam: Optional structure containing the parameters for ODF
%             reconstruction. See the help of reconCSAODF.
% 
% optPerm:    Optimal permutation.
% optFlip:    Optimal flipping.
% grad (out): Corrected gradient table.
% allErrors:  Fiber-continuity error for each permutation/flip configuration.
%
% See also:   reconCSAODF, showODFs, CSAODF_CLI, EXAMPLE, EXAMPLE_CLI.

% Codes by Iman Aganj.
%
% Reference:
% I. Aganj, "Automatic verification of the gradient table in
% diffusion-weighted MRI based on fiber continuity,"
% Scientific Reports, vol. 8, Article no. 16541, 2018.

function [optPerm, optFlip, grad, allErrors] = checkGradTable(S, S0, grad, Mask, indexGPU, voxelSize, res, reconParam)

if ~exist('res', 'var') || isempty(res)
    res = .5;
end
if ~exist('voxelSize', 'var') || isempty(voxelSize)
    voxelSize = [1 1 1];
end
reconParam.indexGPU = 0;
if ~exist('indexGPU', 'var')
    indexGPU = 0;
end
nDir = 0;
for theta = res:res:(pi-res)
    for phi = 0:(res/sin(theta)):2*pi
        nDir = nDir + 1;
        sampleAngles(nDir,:) = [theta phi 1];
    end
end
S = double(S);
if ~isempty(S0)
    S = S ./ double(S0);
end
if indexGPU
    gpu = gpuDevice(indexGPU);
    disp(['Using the "' gpu.Name '" GPU...']);
    S = gpuArray(S);
    grad = gpuArray(grad);
    sampleAngles = gpuArray(sampleAngles);
end
if ~exist('Mask', 'var') || isempty(Mask)
    warning('No mask was specified. The results may not be precise!')
    if ~isfield(reconParam, 'bValue')
        error('If Mask is not given, then the b-value must be specified in reconParam.bValue.')
    end
    reconParamMask = reconParam;
    reconParamMask.basisOrder = 2;
    Mask = mean(-log(S)/reconParam.bValue, 4)<.001 & makeGFA(reconCSAODF(S, [], cart2sph_phys(grad), reconParamMask))>.4;
elseif indexGPU
    Mask = gpuArray(Mask>.5);
end
Perms = perms(1:3);
negMat = [1 1 1; ones(3) - 2*eye(3)];
s = size(S); s = s(1:3);
Mask = Mask(1:end-1, 1:end-1, 1:end-1);
Mask = Mask(:);
disp('Reconstructing the ODFs.')
sh = reconCSAODF(S, [], cart2sph_phys(grad), reconParam);
clear S
ODF = sampleODFs(sh, sampleAngles);
clear sh
diffODF = zeros([nnz(Mask), nDir, 3], 'like', ODF);
for d = 1:3
    diffODF0 = diff(ODF,1,d) / voxelSize(d);
    diffODF0 = reshape(diffODF0(1:s(1)-1, 1:s(2)-1, 1:s(3)-1, :), [prod(s-1) nDir]);
    diffODF(:,:,d) = diffODF0(Mask, :);
end
allErrors = zeros(size(Perms,1), size(negMat,1), 'like', ODF);
clear ODF diffODF0 Mask S0
disp(['Sampling in ' num2str(nDir) ' directions.'])
diffDirect0 = cart2sph_phys(sampleAngles, true);
for nPerm = 1:size(Perms,1)
    for nNeg = 1:size(negMat,1)
        diffDirect = permute(diffDirect0(:,Perms(nPerm,:)) .* negMat(nNeg,:), [3 1 2]);
        for nDiffDirect = 1:nDir
            allErrors(nPerm, nNeg) = allErrors(nPerm, nNeg) + mean(sum(diffODF(:,nDiffDirect,:) .* diffDirect(1,nDiffDirect,:), 3).^2);
        end
        allErrors(nPerm, nNeg) = allErrors(nPerm, nNeg) / nDir;
        disp(['For permutation [' num2str(Perms(nPerm,:)) '] and direction [' num2str(negMat(nNeg,:)) '], the error is ' num2str(allErrors(nPerm, nNeg)) ' .'])
    end
end
[~, ind] = min(allErrors(:));
[nPerm, nNeg] = ind2sub(size(allErrors), ind);
optPerm = Perms(nPerm,:);
optFlip = negMat(nNeg,:);
if isequal(optPerm, 1:3) && isequal(optFlip, [1 1 1])
    disp('No need to change the gradient table.')
else
    disp('Correct gradient table:')
    disp(['grad = grad(:, [' num2str(optPerm) ']) .* [' num2str(optFlip) '];'])
end
grad = gather(grad(:,optPerm) .* optFlip);
allErrors = gather(allErrors);
