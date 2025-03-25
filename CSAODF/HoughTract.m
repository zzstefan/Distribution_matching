% Performs the Hough-transform based tractography from diffusion MRI ODFs.
%
% tracts = HoughTract(sh, pVol, seedPoints, box2vox, param)
%  
% sh:          4D array of ODFs in the spherical harmonic basis, which can
%              be computed using reconCSAODF.m or reconCSAODF3Q.m.
% pVol:        3D probability map of the existence of fiber bundles. For
%              instance, it can be the GFA (that may be created using
%              makeGFA.m) or FA inside of the white+gray matter mask, and 0
%              outside of it.
% seedPoints:  Seed points that are generated using genSeedPoints.m.
% box2vox:     Optional 3x4 affine transformation matrix that takes every
%              point inside the tractography box [-1,1]^3 to the image
%              voxel space, and can be generated using makeBox2vox.m. If it
%              is [] or not specified, it will be created assuming the
%              image to have isotropic voxels.
% param:       Optional structure containing the following parameters. A
%              default value is set to any non-existing field.
%    param.order:     Number of polynomial coefficients to be optimized for
%                     each angle representing the tangent to the curve
%                     (default: 3).
%    param.lambda:    Prior on length; a larger value results in longer
%                     tracts. If not specified, a suggested value will be
%                     obtained using suggestLambda.m.
%    param.resH:      Discretization resolution (step size) of the Hough
%                     transform parameter space (default: 0.4).
%    param.resS:      Discretization resolution of the tract arc length in
%                     the [-1,1]^3 box (default: 0.01).
%    param.maxS:      Maximum partial tract length (in each side of the
%                     seed point) in the [-1,1]^3 box (default: 1).
%    param.nCores:    Number of CPU cores to be used (default: 1).
%                     Choose Inf to use all the available cores.
%    param.indexGPU:  If provided, the GPU device specified by this index
%                     will be used (default: 0 for no GPU).
%                     Run 'CSAODF_CLI GPUs' to see the indices of the
%                     available GPUs. If N>1 GPUs are available, the user
%                     can specify a 1xN vector of GPU indices for
%                     'param.indexGPU' in addition to param.nCores=N, so
%                     each CPU core uses a different GPU.
%    param.nParts:    If provided, the Hough transform parameter space will
%                     be divided into 'param.nParts' parts, and the code
%                     will only search the part number indicated by
%                     'param.nPart', which needs to be mandatorily
%                     specified (default: 1). This is used by 'CSAODF_CLI
%                     HoughTract' with the 'createJobs' option to create
%                     batch jobs for distributed computing.
%    param.nPart:     See the description for 'param.nParts' above (default: 1).
%    param.verbose:   If more than 0, an estimate of the remaining time
%                     will be displayed every 'param.verbose' iterations
%                     (default: 1000).
% 
% tracts:      A structure containing the tractography results in the
%              following fields:
%    tracts.polyCoeffs:      Optimal polynomial coefficients for each seed point.
%    tracts.partialLengths:  Optimal partial lengths for each seed point.
%    tracts.score:           Score of the chosen tract for each seed point.
%    tracts.param:           A structure containing the above input
%                            parameters, with the additional fields of 'sz'
%                            (size of the ODF array), and 'iter1' and
%                            'iter2' (the start and end iteration numbers
%                            corresponding to the particular 'param.nPart').
%    tracts.seedPoints:      Input seed points.
%    tracts.box2vox:         Input (or automatically created) box2vox.
%
% See also:  reconCSAODF, reconCSAODF3Q, makeGFA, genSeedPoints,
%            makeBox2vox, suggestLambda, showTracts, tracts2vox, connMatrix,
%            exportTrackVis, exportText, sampleTracts, fiberDensityVol,
%            CSAODF_CLI, EXAMPLE, EXAMPLE_CLI.

% Codes by Iman Aganj.
%
% Reference:
% I. Aganj, C. Lenglet, N. Jahanshad, E. Yacoub, N. Harel, P. Thompson, and
% G. Sapiro, "A Hough transform global probabilistic approach to
% multiple-subject diffusion MRI tractography,"
% Medical Image Analysis, vol. 15, no. 4, pp. 414-425, 2011.

function tracts = HoughTract(sh, pVol, seedPoints, box2vox, param)

if min(pVol(:))<0
    error('pVol cannot be negative.')
end
if ~exist('box2vox', 'var') || isempty(box2vox)
    box2vox = makeBox2vox(size(pVol), ones(1,3));
end
if ~exist('param', 'var')
    param = [];
end
if ~isfield(param, 'order')
    param.order = 3;
end
if ~isfield(param, 'resH')
    param.resH = .4;
end
if ~isfield(param, 'resS')
    param.resS = .01;
end
if ~isfield(param, 'maxS')
    param.maxS = 1;
end
if ~isfield(param, 'lambda')
    param.lambda = suggestLambda(sh, pVol);
    disp(['Choosing lambda = ' num2str(param.lambda) ' ...'])
end
if ~isfield(param, 'indexGPU')
    param.indexGPU = 0;
end
if ~isfield(param, 'nParts')
    param.nParts = 1;
    param.nPart = 1;
end
if ~isfield(param, 'verbose')
    param.verbose = 1000;
end
oldVer = verLessThan('matlab', '9.1') || sum(param.indexGPU); % Matlab releases earlier than R2016b.
sz = size(sh);
param.sz = sz;
[tracts.param, tracts.seedPoints, tracts.box2vox] = deal(param, seedPoints, box2vox);
[order, resH, resS, maxS, lambda, indexGPU, nParts, nPart, verbose] = deal(param.order, param.resH, param.resS, param.maxS, param.lambda, param.indexGPU, param.nParts, param.nPart, param.verbose);
if isfield(param, 'nCores')
    nCores = param.nCores;
    if isinf(nCores)
        cl = parcluster;
        nCores = cl.NumWorkers;
    end
    if nCores>1
        disp(['Using ' num2str(nCores) ' cores.'])
        pool = gcp('nocreate');
        if isempty(pool)
            parpool(nCores);
        elseif pool.NumWorkers < nCores
            pool.delete
            parpool(nCores);
        end
        parfor k = 1:nCores
            param2 = rmfield(param , 'nCores');
            param2.nParts = nParts * nCores;
            param2.nPart = (nPart-1)*nCores + k;
            param2.indexGPU = indexGPU(mod(k-1,length(indexGPU))+1); % get(getCurrentTask(), 'ID'))
            partialTracts(k) = HoughTract(sh, pVol, seedPoints, box2vox, param2);
        end
        sz2 = size(partialTracts(1).polyCoeffs);
        polyCoeffsAll = zeros([sz2 nCores]);
        partialLengthsAll = zeros(sz2(1),2,nCores);
        polyCoeffsAll(:) = [partialTracts.polyCoeffs];
        partialLengthsAll(:) = [partialTracts.partialLengths];
        [score, ind] = max([partialTracts.score],[],2);
        polyCoeffs = zeros(sz2(1:2));
        partialLengths = zeros(sz2(1),2);
        for iX = 1:sz2(1)
            polyCoeffs(iX,:) = polyCoeffsAll(iX,:,ind(iX));
            partialLengths(iX,:) = partialLengthsAll(iX,:,ind(iX));
        end
        [tracts.partialLengths, tracts.polyCoeffs, tracts.score] = deal(partialLengths, polyCoeffs, score);
        return
    end
end

mask = pVol>0;
if indexGPU
    gpu = gpuDevice(indexGPU);
    dispStr = ['Using the "' gpu.Name '" GPU...'];
    if nParts>1
        dispStr = [dispStr ' (part ' num2str(nPart) '/' num2str(nParts) ')'];
    end
    disp(dispStr)
    seedPoints = gpuArray(seedPoints);
    sh = gpuArray(sh);
    pVol = gpuArray(pVol);
    mask = gpuArray(mask);
    box2vox = gpuArray(box2vox);
end
seqS = 0:resS:maxS;
nSeqS = length(seqS);
xCurve = zeros(nSeqS, 3, 1, 2);
baseSeqH = 0:resH:1.6;
baseSeqH = [-baseSeqH(end:-1:2) baseSeqH];
nSeqH = length(baseSeqH);
nTotalIters = nSeqH^(2*order);
nSeqHones = nSeqH*ones(1,2*order);
siSize = [2*order nSeqH];
kModN = [1:order, 1:order]';
seqH = ((2-1./kModN)./(maxS.^(kModN-1))) * baseSeqH;
basisOrder = (sqrt(8*sz(end)+1)-3)/2;
nSeedPoints = size(seedPoints, 1);
partialLengths = zeros(nSeedPoints, 2, 'like', sh);
score = -inf(nSeedPoints, 1, 'like', sh);
H0 = -inf(nSeqS, nSeedPoints, 2, 'like', sh);
indS = int32(ndgrid(1:nSeqS, 1:nSeedPoints));
indS(:,:,2) = indS(:,:,1) + nSeqS;
sz = int32(sz);
sub2indCoeff = [sz(1), sz(1)*sz(2)];
if indexGPU
    indS = gpuArray(indS);
    sub2indCoeff = gpuArray(sub2indCoeff);
end
polyCoeffCell = cell(1,2*order);

seedPoints = permute(seedPoints, [3 2 1]);
sh = reshape(sh, [], sz(end));
if oldVer
    sh = bsxfun(@times, sh, precompSH(basisOrder)');
    for iDir = 1:2
        sPower(:,:,iDir) = bsxfun(@power, (2*iDir-3)*seqS, (0:order-1)');
    end
else
    sh = sh .* precompSH(basisOrder)';
    for iDir = 1:2
        sPower(:,:,iDir) = ((2*iDir-3)*seqS) .^ ((0:order-1)');
    end
end
pVol = log(double(pVol)) + lambda;
indH = ones(nSeedPoints,1);
iter1 = floor(nTotalIters*(nPart-1)/nParts)+1;
iter2 = min(floor(nTotalIters*nPart/nParts)-1, nTotalIters)+1;
param.iter1 = iter1;
param.iter2 = iter2;
for i = iter1:iter2
    if verbose && mod(i,verbose)==0
        dispStr = ['Iteration ' num2str(i) ' out of ' num2str(nTotalIters) '.'];
        if nParts>1
            dispStr = [dispStr ' (part ' num2str(nPart) '/' num2str(nParts) ')'];
        end
        if exist('tStart', 'var')
            t = toc(tStart);
            dispStr = sprintf('%s Remaining: ~ %.1f min.', dispStr, t*(iter2-i)/(60*verbose));
        end
        disp(dispStr)
        tStart = tic;
    end
    [polyCoeffCell{:}] = ind2sub(nSeqHones, i);
    polyCoeff = cell2mat(polyCoeffCell);
    h = seqH(sub2ind(siSize, 1:(2*order), polyCoeff));
    Y = zeros(2*nSeqS, sz(end));
    for iDir = 1:2
        theta = h(1:order) * sPower(:, :, iDir);
        phi = h(order+1:2*order) * sPower(:, :, iDir);
        sTheta = sin(theta);
        cTheta = cos(theta);
        Y((1:nSeqS)+(iDir-1)*nSeqS,:) = compSH(basisOrder, theta, phi, [], cTheta, sTheta)';
        xCurve(:,:,1,iDir) = [0 0 0; ((2*iDir-3)*resS)*[sTheta(1:end-1).*cos(phi(1:end-1)); sTheta(1:end-1).*sin(phi(1:end-1)); cTheta(1:end-1)]'];
    end
    if indexGPU
        Y = gpuArray(Y);
    end
    if oldVer
        x = bsxfun(@plus, cumsum(xCurve), seedPoints);
        X1 = int32(sum(bsxfun(@times, x, box2vox(1,1:3)),2) + box2vox(1,4));
        X2 = int32(sum(bsxfun(@times, x, box2vox(2,1:3)),2) + box2vox(2,4));
        X3 = int32(sum(bsxfun(@times, x, box2vox(3,1:3)),2) + box2vox(3,4));
    else
        x = cumsum(xCurve) + seedPoints;
        X1 = int32(sum(x .* box2vox(1,1:3),2) + box2vox(1,4));
        X2 = int32(sum(x .* box2vox(2,1:3),2) + box2vox(2,4));
        X3 = int32(sum(x .* box2vox(3,1:3),2) + box2vox(3,4));
    end
    legitX = X1>0 & X2>0 & X3>0 & X1<=sz(1) & X2<=sz(2) & X3<=sz(3);
    indX = X1(legitX) + sub2indCoeff(1)*X2(legitX) + sub2indCoeff(2)*X3(legitX) -(sub2indCoeff(1)+sub2indCoeff(2));
    legitX0 = legitX;
    legitX(legitX) = mask(indX);
    legitX = cummin(legitX);
    indX = indX(legitX(legitX0));
    odf = sum(sh(indX,:) .* Y(indS(legitX),:), 2);
    H = H0;
    H(legitX) = log(max(odf, 0) + .02) + pVol(indX);
    [maxPartialH, indPartialH] = max(cumsum(H));
    maxH = sum(maxPartialH,3)';
    xUpdate = maxH > score;
    partialLengths(xUpdate,:) = indPartialH(1,xUpdate,:);
    score(xUpdate) = maxH(xUpdate);
    indH(xUpdate) = i;
end
[polyCoeffCell{:}] = ind2sub(nSeqHones, indH);
polyCoeff = cell2mat(polyCoeffCell);
polyCoeffs  = zeros(nSeedPoints, 2*order);
for iPoint = 1:nSeedPoints
    polyCoeffs(iPoint,:) = seqH(sub2ind(siSize, 1:(2*order), polyCoeff(iPoint,:)));
end
if indexGPU
    partialLengths = gather(partialLengths);
    score = gather(score);
end
[tracts.partialLengths, tracts.polyCoeffs, tracts.score, tracts.param] = deal(partialLengths, polyCoeffs, score, param);
