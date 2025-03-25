% Suggests an optimal value for the length prior 'lambda', based on the
% histogram of the local score that is computed using the 4D array 'sh' of
% ODFs in the spherical harmonic basis and the 3D probability map 'pVol' of
% the existence of fiber bundles. If the optional 'showHist' is true, the
% histogram is plotted.
%
% lambda = suggestLambda(sh, pVol, showHist)
%  
% See also:  HoughTract, genSeedPoints, reconCSAODF, reconCSAODF3Q,
%            sampleODFs, CSAODF_CLI, EXAMPLE, EXAMPLE_CLI.

% Code by Iman Aganj.

function lambda = suggestLambda(sh, pVol, showHist)

angles = cart2sph_phys(eye(3));
sh = reshape(sh, [], 1, 1, size(sh,4));
sh = sh(pVol(:)>0,1,1,:);
pVol = pVol(pVol(:)>0);
odf = sampleODFs(sh, angles);
score = bsxfun(@plus, log(odf.*(odf>0) + .02), log(pVol));
score = score(:);

if exist('showHist', 'var') && showHist
    figure
    hist(score, 100)
    xlabel(['min: ' num2str(min(score)) ',    max: ' num2str(max(score))])
end
score = score(score>min(score));
lambda = -(mean(score)+std(score));
