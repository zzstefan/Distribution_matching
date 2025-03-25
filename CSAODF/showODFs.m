% Visualizes the ODFs. An optional parameter will be ignored if [] is
% entered for it.
%
% showODFs(sh, ROI, overlayVol, Mask, gridSize, scaleODF, dnSample, showNegative, directColor, figPosition)
%  
% sh:            4D array of ODFs in the spherical harmonic basis.
% ROI:           3x2 matrix [x1,x2;y1,y2;z1,z2] that specifies the corners
%                of the region of interest. It is best to choose the size
%                of the cuboid to be 1 in at least one of the dimensions.
% overlayVol:    Optional 3D image over which the ODFs will be visualized
%                (default: GFA). Enter 0 for no overlay.
% Mask:          If provided, voxels where Mask is 0 will not be shown.
% gridSize:      Optional resolution of the ODF surface (default: 30).
% scaleODF:      Optional size of the visualized ODFs (default: 1.5).
% dnSample:      Optional downsampling rate to visualize faster (default: 1).
% showNegative:  If true, negative values will be shown (as black) (default: false).
% directColor:   If true, directional color-coding will be used (default: false).
% figPosition:   Optional figure 'Position' property.
%
% See also:  reconCSAODF, reconCSAODF3Q, sampleODFs, makeGFA, showTracts,
%            CSAODF_CLI, EXAMPLE, EXAMPLE_CLI.

% Code by Iman Aganj.

function fig = showODFs(sh, ROI, overlayVol, Mask, gridSize, scaleODF, dnSample, showNegative, directColor, figPosition)

if ~exist('overlayVol', 'var') || isempty(overlayVol)
    overlayVol = makeGFA(sh);
end
sz = size(sh);
if ~exist('Mask', 'var') || isempty(Mask)
    Mask = true(sz(1:3));
end
if ~exist('gridSize', 'var') || isempty(gridSize)
    gridSize = 30;
end
if ~exist('scaleODF', 'var') || isempty(scaleODF)
    scaleODF = 1.5;
end
if ~exist('dnSample', 'var') || isempty(scaleODF)
    dnSample = 1;
end
if ~exist('directColor', 'var') || isempty(directColor)
    directColor = false;
end
oldVer = verLessThan('matlab', '9.1'); % Matlab releases earlier than R2016b.
basisOrder = (sqrt(8*sz(end)+1)-3)/2;
[x,y,z] = sphere(gridSize);
angles = cart2sph_phys([x(:) y(:) z(:)]);
sTheta = sin(reshape(angles(:,1),size(x)));
[Y, dY] = compSH(basisOrder, angles(:,1)', angles(:,2)');
Y = permute(Y, [3 4 5 1 2]);
dY = permute(dY, [4 5 6 1 2 3]);
sh = sh(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2),ROI(3,1):ROI(3,2),:);
Mask = Mask(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2),ROI(3,1):ROI(3,2));
if dnSample > 1
    MaskDnSample = false(size(Mask));
    MaskDnSample(1:dnSample:end, 1:dnSample:end, 1:dnSample:end) = true;
    Mask = Mask & MaskDnSample;
    clear MaskDnSample
    scaleODF = scaleODF * dnSample;
end
sh = sh * scaleODF;
if oldVer
    odf = sum(bsxfun(@times, sh, Y), 4);
else
    odf = sum(sh .* Y, 4);
end
if ~(exist('showNegative', 'var') && ~isempty(showNegative) && showNegative)
    odf = odf .* (odf>0);
end
odf = permute(reshape(odf, [size(sh,1) size(sh,2) size(sh,3) size(x)]), [4 5 1 2 3]);
if oldVer
    odfX = bsxfun(@times, odf, cat(6, x, y, z));
    dOdf = sum(bsxfun(@times, sh, dY), 4);
else
    odfX = odf .* cat(6, x, y, z);
    dOdf = sum(sh .* dY, 4);
end
dOdf = permute(reshape(dOdf, [size(sh,1) size(sh,2) size(sh,3) size(x) 2]), [4 5 6 1 2 3]);
rH = cat(3, x, y, z);
rho = sqrt(x.^2+y.^2);
thH = cat(3, x.*z./rho, y.*z./rho, -rho);
phH = cat(3, -y./rho, x./rho, zeros(size(x)));
cMap = jet(1000);
odfVec = zeros(numel(x), 1);
if exist('figPosition', 'var') && ~isempty(figPosition)
    fig = figure('Position', figPosition);
else
    fig = figure;
end
hold on
axis equal %tight
set(gca, 'Clipping', 'off')
indSing = find(ROI(:,1)==ROI(:,2),1);
if ~isempty(indSing)
    switch indSing
        case 1
            view(90,0)
        case 2
            view(0,0)
    end
end
nODFs = nnz(Mask);
nODF = 0;
nODFdraw = 0;
nTotalDraw = 10;
for i = 1:size(odf,3)
    for j = 1:size(odf,4)
        for k = 1:size(odf,5)
            if Mask(i,j,k)
                nODF = nODF + 1;
                odfVec(:) = odf(:,:,i,j,k);
                indPos = odfVec>=0;
                odfVec = odfVec - min(odfVec(indPos));
                if directColor
                    cData = zeros(length(indPos),3);
                    if oldVer
                        odfXVec = reshape(bsxfun(@rdivide, abs(odfX(:,:,i,j,k,:)), sqrt(sum(odfX(:,:,i,j,k,:).^2,6))), size(cData));
                    else
                        odfXVec = reshape(abs(odfX(:,:,i,j,k,:)) ./ sqrt(sum(odfX(:,:,i,j,k,:).^2,6)), size(cData));
                    end
                    cData(indPos,:) = odfXVec(indPos,:);
                    cData = reshape(cData, size(x,1), size(x,2), 3);
                else
                    if any(odfVec>0)
                        cData = zeros(length(indPos),3);
                        odfVec = odfVec / max(odfVec);
                        cData(indPos,:) = cMap(floor(odfVec(indPos)*(size(cMap,1)-1))+1, :);
                        cData = reshape(cData, size(x,1), size(x,2), 3);
                    else
                        cData = repmat(permute(cMap(round((size(cMap,1)+1)/2),:), [3 1 2]), size(x,1), size(x,2));
                    end
                end
                hS = surf(odfX(:,:,i,j,k,1)+ROI(1,1)+i-1, odfX(:,:,i,j,k,2)+ROI(2,1)+j-1, odfX(:,:,i,j,k,3)+ROI(3,1)+k-1, cData, 'EdgeColor', 'none', 'FaceColor', 'interp');
                if oldVer
                    vertexNormals = rH - bsxfun(@times, dOdf(:,:,1,i,j,k)./odf(:,:,i,j,k), thH) - bsxfun(@times, dOdf(:,:,2,i,j,k)./(odf(:,:,i,j,k).*sTheta), phH);
                    vertexNormals = bsxfun(@rdivide, vertexNormals, sqrt(sum(vertexNormals.^2,3)));
                else
                    vertexNormals = rH - (dOdf(:,:,1,i,j,k)./odf(:,:,i,j,k)) .* thH - (dOdf(:,:,2,i,j,k)./(odf(:,:,i,j,k).*sTheta)) .* phH;
                    vertexNormals = vertexNormals ./ sqrt(sum(vertexNormals.^2,3));
                end
                hS.VertexNormals = vertexNormals;
                hS.FaceLighting = 'gouraud';
                hS.AmbientStrength = .8;
                if round(nTotalDraw*nODF/nODFs) > nODFdraw
                    set(fig, 'Name', ['(' num2str(round(100*nODF/nODFs)) '%)'])
                    drawnow
                    nODFdraw = round(nTotalDraw*nODF/nODFs);
                end
            end
        end
    end
end
set(fig, 'Name', '')
camlight
if numel(overlayVol)>1 || overlayVol
    overlayVol = permute(overlayVol, [2 1 3]);
    sl = mean(ROI,2);
    for d = 1:3
        if size(overlayVol, d) == 1
            overlayVol = cat(d, overlayVol, overlayVol);
        end
    end
    try
        hSlice = slice(overlayVol, sl(1), sl(2), sl(3));
        for i = 1:3
            V = hSlice(i).CData;
            V = V - min(V(:)); V = V / max(V(:));
            hSlice(i).CData = repmat(V, [1 1 3]);
            hSlice(i).FaceColor = 'interp';
            hSlice(i).FaceAlpha = .5;
            hSlice(i).AmbientStrength = 1;
        end
    catch ME
        warning(ME.message)
    end
end
shading interp
axis equal
xlabel X, ylabel Y, zlabel Z
hold off
drawnow
if nargout==0
    clear fig
end
