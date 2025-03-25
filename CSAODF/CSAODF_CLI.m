% This is the command-line interface for CSA-ODF reconstruction,
% Hough-transform tractography, and visualization of the results. It is
% especially useful to compile this file with 'mcc -m CSAODF_CLI', so it
% can be run on computers with no Matlab license. Run CSAODF_CLI without
% arguments to see the help page.
%
% See also:   EXAMPLE_CLI, EXAMPLE, reconCSAODF, HoughTract.

% Codes by Iman Aganj

function CSAODF_CLI(varargin)

if nargin==0
    varargin{1} = 'none';
end
p = inputParser;
p.FunctionName = mfilename;
p.StructExpand = false;
switch upper(varargin{1})
    case 'RECONCSAODF'
        if nargin<5
            disp('Computes the q-ball CSA-ODFs from either a single shell or three shells')
            disp('of diffusion MRI data.')
            disp(' ')
            disp('CSAODF_CLI reconCSAODF <input_DWIs> <input_bvecs> <input_bvals>')
            disp('           <output_ODFs> [bValue <bValue>] [basisOrder <basisOrder>]')
            disp('           [delta <delta>] [method <1|2|3>] [roundBVal <roundBVal>]')
            disp('           [indexGPU <indexGPU>]')
            disp(' ')
            disp('input_DWIs:    4D NIFTI file containing the diffusion-weighted images.')
            disp('input_bvecs:   Text file containing the gradient table.')
            disp('input_bvals:   Text file containing the b-values.')
            disp('output_ODFs:   Name of the NIFTI file where the computed ODFs in the')
            disp('               spherical harmonic basis will be saved.')
            disp(' ')
            disp('Optional parameters:')
            disp('bValue:        For single-shell reconstruction, bValue is the b-value of')
            disp('               the desired shell to be used. If not specified, the')
            disp('               smallest b-value found in ''input_bvals'' (after b=0) is')
            disp('               chosen.')
            disp('               For three-shell reconstruction, use')
            disp('               [bValue1,bValue2,bValue3]. Note that the b-values of the')
            disp('               shells must be such that 0,bValue1,bValue2,bValue3 is an')
            disp('               arithmetic sequence, and the gradient orientations must be')
            disp('               the same for all three shells.')
            disp('basisOrder:    Spherical harmonic basis order (an even number, default: 4).')
            disp('delta:         Regularization parameter (default: 0.0001).')
            disp('method:        Reconstruction method:')
            disp('                  1 (default) --> CSA-ODF')
            disp('                  2           --> Original (Tuch''s) ODF')
            disp('                  3           --> Original (Tuch''s) ODF, with')
            disp('                                  Laplace-Beltrami sharpening (lambda=0.1)')
            disp('roundBVal:     Rounds the b-values in ''input_bvals''. Similar to the')
            disp('               ''round'' command in Matlab, if ''roundBVal'' is negative, the')
            disp('               b-values are rounded abs(''roundBVal'') digits to the left')
            disp('               of the decimal point. For example, if ''roundBVal'' is -2,')
            disp('               the b-value of 4990 will be rounded to 5000.')
            disp('indexGPU:      If provided, the GPU device specified by this index will')
            disp('               be used. Run ''CSAODF_CLI GPUs'' to see the indices of the')
            disp('               available GPUs.')
            disp(' ')
            disp('Reference:')
            disp('I. Aganj, C. Lenglet, G. Sapiro, E. Yacoub, K. Ugurbil, and N. Harel,')
            disp('"Reconstruction of the orientation distribution function in single and')
            disp('multiple shell q-ball imaging within constant solid angle,"')
            disp('Magnetic Resonance in Medicine, vol. 64, no. 2, pp. 554-566, 2010.')
            return
        end
        p.addRequired('input_DWIs');
        p.addRequired('input_bvecs');
        p.addRequired('input_bvals');
        p.addRequired('output_ODFs');
        p.addParameter('bValue', []);
        p.addParameter('basisOrder', []);
        p.addParameter('delta', []);
        p.addParameter('method', []);
        p.addParameter('indexGPU', []);
        p.addParameter('roundBVal', []);
        p.parse(varargin{2:end});
        bvecs = load(p.Results.input_bvecs, 'ascii');
        bvals = load(p.Results.input_bvals, 'ascii');
        if size(bvecs,2)>size(bvecs,1)
            bvecs = bvecs';
        end
        if size(bvals,2)>size(bvals,1)
            bvals = bvals';
        end
        if ~isempty(p.Results.roundBVal)
            bvals = round(bvals, str2double(p.Results.roundBVal));
        end
        b = unique(bvals);
        disp(array2table([b sum(bsxfun(@eq, b(:,1), bvals'),2)], 'VariableNames', {'bValue', 'numberOfImages'}))
        if isempty(p.Results.bValue)
            bValue = b(2);
        else
            bValue = str2num(p.Results.bValue);
        end
        nShells = length(bValue);
        bValue = sort(bValue, 'ascend');
        bValueExt = [0 bValue];
        if nShells==3
            if any(diff(bValueExt, 2))
                warning(['The series 0,bValue(1),bValue(2),bValue(3) must be an arithmetic sequence. 0,' num2str(bValue(1)) ',' num2str(bValue(2)) ',' num2str(bValue(3)) ' is not. The results may not be accurate!'])
            end
        elseif nShells~=1
            error('bValue must contain either 1 element or 3 elements, respectively for 1-shell and 3-shell reconstructions.')
        end
        indB = arrayfun(@(x) find(bvals==x), bValueExt', 'UniformOutput', false);
        for i = 1:nShells+1
            if isempty(indB{i})
                error(['No b=' num2str(bValueExt(i)) ' images were found.'])
            end
            bvecsShell{i} = bvecs(indB{i},:);
            if i>2 && any(bvecsShell{i}(:) ~= bvecsShell{i-1}(:))
                warning(['The gradient tables of b=' num2str(bValueExt(i-1)) ' and b=' num2str(bValueExt(i)) ' are not identical. The results may not be accurate!'])
            end
        end
        param = [];
        if ~isempty(p.Results.basisOrder)
            param.basisOrder = str2double(p.Results.basisOrder);
        end
        if ~isempty(p.Results.delta)
            param.delta = str2double(p.Results.delta);
        end
        if ~isempty(p.Results.method)
            if nShells==1
                param.method = str2double(p.Results.method);
            else
                warning('''Method'' is ignored in three-shell reconstruction.')
            end
        end
        if ~isempty(p.Results.indexGPU)
            param.indexGPU = str2double(p.Results.indexGPU);
        end
        angles = cart2sph_phys(bvecsShell{2});
        disp('Reading the images...')
        DWIs = load_nifti_allOS(p.Results.input_DWIs);
        DWIs.vol = DWIs.vol(:,:,:,cell2mat(indB));
        newIndB(1,:) = [1 length(indB{1})];
        for i = 2:nShells+1
            newIndB(i,:) = newIndB(i-1,2) + [1 length(indB{i})];
        end
        disp(['Using b=' num2str(bValue, '%d,') ' for ' num2str(nShells) '-shell ODF reconstruction...'])
        if nShells==1
            DWIs.vol = reconCSAODF(DWIs.vol(:,:,:,newIndB(2,1):newIndB(2,2)), mean(DWIs.vol(:,:,:,newIndB(1,1):newIndB(1,2)),4), angles, param);
        else
            DWIs.vol = reconCSAODF3Q(DWIs.vol(:,:,:,newIndB(2,1):newIndB(2,2)), DWIs.vol(:,:,:,newIndB(3,1):newIndB(3,2)), DWIs.vol(:,:,:,newIndB(4,1):newIndB(4,2)), mean(DWIs.vol(:,:,:,newIndB(1,1):newIndB(1,2)),4), angles, param);
        end
        DWIs.dim(5) = size(DWIs.vol,4);
        DWIs.datatype = 64;
        save_nifti_allOS(DWIs, p.Results.output_ODFs);
    case 'SHOWODFS'
        if nargin<2
            disp('Visualizes the ODFs.')
            disp(' ')
            disp('CSAODF_CLI showODFs <input_ODFs> [overlay <overlay>] [Mask <Mask>]')
            disp('           [ROI <ROI>] [gridSize <gridSize>] [scaleODF <scaleODF>]')
            disp('           [dnSample <dnSample>] [showNegative <showNegative>]')
            disp('           [directColor <directColor>]')
            disp(' ')
            disp('input_ODFs:     4D NIFTI file with ODFs in the spherical harmonic basis.')
            disp(' ')
            disp('Optional parameters:')
            disp('overlay:        Image (in NIFTI), over which the ODFs will be visualized')
            disp('                (default: GFA). ')
            disp('Mask:           Mask image (in NIFTI); voxels where Mask is 0 will not be shown.')
            disp('ROI:            Corners of the region of interest, [x1,x2,y1,y2,z1,z2].')
            disp('                It is best to choose the size of the cuboid to be 1 in at')
            disp('                least one of the dimensions.')
            disp('gridSize:       Resolution of the ODF surface (default: 30).')
            disp('scaleODF:       Size of the visualized ODFs (default: 1.5).')
            disp('dnSample:       Optional downsampling rate to visualize faster (default: 1).')
            disp('showNegative:   If 1, negative values will be shown (as black) (default: 0).')
            disp('directColor:    If 1, directional color-coding will be used (default: 0).')
            return
        end
        p.addRequired('input_ODFs');
        p.addParameter('overlay', []);
        p.addParameter('Mask', []);
        p.addParameter('ROI', []);
        p.addParameter('gridSize', []);
        p.addParameter('scaleODF', []);
        p.addParameter('dnSample', []);
        p.addParameter('showNegative', []);
        p.addParameter('directColor', '0');
        p.parse(varargin{2:end});
        sh = load_nifti_allOS(p.Results.input_ODFs);
        if isempty(p.Results.overlay)
            overlayVol = [];
        else
            overlayVol = load_nifti_allOS(p.Results.overlay);
            overlayVol = overlayVol.vol;
        end
        if isempty(p.Results.Mask)
            Mask = [];
        else
            Mask = load_nifti_allOS(p.Results.Mask);
            Mask = Mask.vol;
        end
        if isempty(p.Results.ROI)
            sz = size(sh.vol);
            ROI = round(bsxfun(@times, sz(1:3)', [.3 .7; .5 .5; .2 .9]));
        else
            ROI = str2num(p.Results.ROI);
            if numel(ROI)~=6
                error('ROI needs to have exactly 6 elements, in the format [x1,x2,y1,y2,z1,z2].')
            end
            if size(ROI,1)==1
                ROI = [ROI(1:2); ROI(3:4); ROI(5:6)];
            end
        end
        if isempty(p.Results.gridSize)
            gridSize = [];
        else
            gridSize = str2double(p.Results.gridSize);
        end
        if isempty(p.Results.scaleODF)
            scaleODF = [];
        else
            scaleODF = str2double(p.Results.scaleODF);
        end
        if isempty(p.Results.dnSample)
            dnSample = [];
        else
            dnSample = str2double(p.Results.dnSample);
        end
        if isempty(p.Results.showNegative)
            showNegative = [];
        else
            showNegative = str2double(p.Results.showNegative);
        end
        showODFs(sh.vol, ROI, overlayVol, Mask, gridSize, scaleODF, dnSample, showNegative, str2double(p.Results.directColor))
    case 'VERIFY'
        if nargin<5
            disp('Verifies if the gradient table uses the same coordinates as the image.')
            disp(' ')
            disp('CSAODF_CLI verify <input_DWIs> <input_bvecs> <input_bvals> <output_bvecs>')
            disp('           [mask <Mask>] [res <res>] [*other reconCSAODF options]')
            disp(' ')
            disp('input_DWIs:    4D NIFTI file containing the diffusion-weighted images.')
            disp('input_bvecs:   Text file containing the gradient table.')
            disp('input_bvals:   Text file containing the b-values.')
            disp('output_bvecs:  Name of the text file where the corrected gradient table')
            disp('               will be saved.')
            disp(' ')
            disp('Optional parameters:')
            disp('mask:          Fibrous tissue mask (in NIFTI). If not specified (not')
            disp('               recommended), a rough mask will be created.')
            disp('res:           Angular resolution in radians (default: 0.5).')
            disp(' ')
            disp('*Others:       Run ''CSAODF_CLI reconCSAODF'' to see other optional ODF-')
            disp('               reconstruction parameters.')
            disp(' ')
            disp('Reference:')
            disp('I. Aganj, "Automatic verification of the gradient table in')
            disp('diffusion-weighted MRI based on fiber continuity,"')
            disp('Scientific Reports, vol. 8, Article no. 16541, 2018.')
            return
        end
        p.addRequired('input_DWIs');
        p.addRequired('input_bvecs');
        p.addRequired('input_bvals');
        p.addRequired('output_bvecs');
        p.addParameter('mask', []);
        p.addParameter('res', []);
        p.addParameter('bValue', []);
        p.addParameter('basisOrder', []);
        p.addParameter('delta', []);
        p.addParameter('method', []);
        p.addParameter('indexGPU', []);
        p.addParameter('roundBVal', []);
        p.parse(varargin{2:end});
        bvecs = load(p.Results.input_bvecs, 'ascii');
        bvals = load(p.Results.input_bvals, 'ascii');
        if size(bvecs,2)>size(bvecs,1)
            bvecs = bvecs';
        end
        if size(bvals,2)>size(bvals,1)
            bvals = bvals';
        end
        if ~isempty(p.Results.roundBVal)
            bvals = round(bvals, str2double(p.Results.roundBVal));
        end
        b = unique(bvals);
        disp(array2table([b sum(bsxfun(@eq, b(:,1), bvals'),2)], 'VariableNames', {'bValue', 'numberOfImages'}))
        if isempty(p.Results.bValue)
            bValue = b(2);
        else
            bValue = str2num(p.Results.bValue);
        end
        if length(bValue)>1
            error('Multi-shell reconstruction is not supported here.')
        end
        bValue = sort(bValue, 'ascend');
        bValueExt = [0 bValue];
        indB = arrayfun(@(x) find(bvals==x), bValueExt', 'UniformOutput', false);
        for i = 1:2
            if isempty(indB{i})
                error(['No b=' num2str(bValueExt(i)) ' images were found.'])
            end
            bvecsShell{i} = bvecs(indB{i},:);
        end
        param = [];
        if ~isempty(p.Results.res)
            res = str2double(p.Results.res);
        else
            res = [];
        end
        if ~isempty(p.Results.mask)
            Mask = load_nifti_allOS(p.Results.mask);
            Mask = Mask.vol > .5;
        else
            Mask = [];
            param.bValue = bValue;
        end
        if ~isempty(p.Results.basisOrder)
            param.basisOrder = str2double(p.Results.basisOrder);
        end
        if ~isempty(p.Results.delta)
            param.delta = str2double(p.Results.delta);
        end
        if ~isempty(p.Results.method)
            param.method = str2double(p.Results.method);
        end
        if ~isempty(p.Results.indexGPU)
            indexGPU = str2double(p.Results.indexGPU);
        else
            indexGPU = 0;
        end
        disp('Reading the images...')
        DWIs = load_nifti_allOS(p.Results.input_DWIs);
        DWIs.vol = DWIs.vol(:,:,:,cell2mat(indB));
        newIndB(1,:) = [1 length(indB{1})];
        newIndB(2,:) = newIndB(1,2) + [1 length(indB{2})];
        disp(['Using b=' num2str(bValue, '%d,') ' for ODF reconstruction...'])
        [optPerm, optFlip] = checkGradTable(DWIs.vol(:,:,:,newIndB(2,1):newIndB(2,2)), mean(DWIs.vol(:,:,:,newIndB(1,1):newIndB(1,2)),4), bvecsShell{2}, Mask, indexGPU, DWIs.pixdim(2:4), res, param);
        dlmwrite(p.Results.output_bvecs, bvecs(:,optPerm).*optFlip, 'delimiter', '\t', 'precision', 8)
    case 'MAKEGFA'
        if nargin<3
            disp('Computes the generalized fractional anisotropy from the ODFs.')
            disp(' ')
            disp('CSAODF_CLI makeGFA <input_ODFs> <output_GFA>')
            disp(' ')
            disp('input_ODFs:   4D NIFTI file with ODFs in the spherical harmonic basis.')
            disp('output_GFA:   Name of the output NIFTI file to save the GFA.')
            return
        end
        p.addRequired('input_ODFs');
        p.addRequired('output_GFA');
        p.parse(varargin{2:end});
        sh = load_nifti_allOS(p.Results.input_ODFs);
        sh.dim(5) = 1;
        sh.vol = makeGFA(sh.vol);
        sh.datatype = 64;
        save_nifti_allOS(sh, p.Results.output_GFA);
    case 'GENSEEDPOINTS'
        if nargin<3
            disp('Generates quasi-random seed points (in the [-1,1]^3 box) for')
            disp('Hough-transform tractography.')
            disp(' ')
            disp('CSAODF_CLI genSeedPoints <output_seedPoints> <seedMask> [pVol <pVol>]')
            disp('           [input_ODFs <input_ODFs>] [numSeedPoints <numSeedPoints>]')
            disp('           [voxelSize <voxelSize>]')
            disp(' ')
            disp('output_seedPoints:   Name of the output .mat file where the seed points')
            disp('                     and the box2vox matrix will be saved. ')
            disp('seedMask:            Brain mask (in NIFTI) that excludes the CSF.')
            disp(' ')
            disp('Optional parameters:')
            disp('pVol:                Probability image (in NIFTI) of the existence of')
            disp('                     fiber bundles. If not provided, the ODFs should be')
            disp('                     given with the ''input_ODFs'' option so the')
            disp('                     generalized fractional anisotropy be computed and')
            disp('                     used as pVol.')
            disp('input_ODFs:          4D NIFTI file containing the ODFs in the spherical')
            disp('                     harmonic basis. It should be given only if ''pVol'' is')
            disp('                     not provided (see the description for ''pVol'' above).')
            disp('numSeedPoints:       Number of seed points to be generated (default: 10000).')
            disp('voxelSize:           1x3 vector, in the format [s1,s2,s3], which contains')
            disp('                     the voxel size. If not provided, it is computed from')
            disp('                     the header of ''seedMask''.')
            return
        end
        p.addRequired('output_seedPoints');
        p.addRequired('seedMask');
        p.addParameter('pVol', []);
        p.addParameter('input_ODFs', []);
        p.addParameter('numSeedPoints', '10000');
        p.addParameter('voxelSize', []);
        p.parse(varargin{2:end});
        seedMask = load_nifti_allOS(p.Results.seedMask);
        if isempty(p.Results.voxelSize)
            vox2ras = [seedMask.srow_x seedMask.srow_y seedMask.srow_z]';
            box2vox = makeBox2vox(size(seedMask.vol), vox2ras);
        else
            box2vox = makeBox2vox(size(seedMask.vol), str2num(p.Results.voxelSize));
        end
        if isempty(p.Results.pVol)
            if isempty(p.Results.input_ODFs)
                error('Either pVol should be specified, or to use the GFA for it, input_ODFs needs to be provided. For uniform seed points, seedMask can be used for pVol.')
            end
            sh = load_nifti_allOS(p.Results.input_ODFs);
            pVol = makeGFA(sh.vol);
            clear sh
        else
            pVol = load_nifti_allOS(p.Results.pVol);
            pVol = pVol.vol;
        end
        seedPoints = genSeedPoints(str2double(p.Results.numSeedPoints), pVol.*(seedMask.vol>0), box2vox);
        save(p.Results.output_seedPoints, 'seedPoints', 'box2vox')
    case 'HOUGHTRACT'
        if nargin<5
            disp('Performs the Hough-transform based tractography from diffusion MRI ODFs.')
            disp(' ')
            disp('CSAODF_CLI HoughTract <input_ODFs> <input_seedPoints> <output_tracts>')
            disp('           <mask> [pVol <pVol>] [lambda <lambda>] [order <order>]')
            disp('           [resH <resH>] [resS <resS>] [maxS <maxS>] [indexGPU <indexGPU>]')
            disp('           [nCores <nCores>] [createJobs <numJobs>]')
            disp('           [jobFilename <jobFilename>] [verbose <numIters>]')
            disp('           [nParts <nParts>] [nPart <nPart>]')
            disp(' ')
            disp('input_ODFs:         4D NIFTI file with ODFs in the spherical harmonic basis.')
            disp('input_seedPoints:   The .mat file containing the seed points and the')
            disp('                    box2vox matrix, generated by ''CSAODF_CLI genSeedPoints''.')
            disp('output_tracts:      Name of the output .mat file where the tractography')
            disp('                    results will be saved. ')
            disp('mask:               Brain mask (in NIFTI) that excludes the CSF.')
            disp(' ')
            disp('Optional parameters:')
            disp('pVol:               Probability image (in NIFTI) of the existence of')
            disp('                    fiber bundles. If not provided, the generalized')
            disp('                    fractional anisotropy will be created and used.')
            disp('lambda:             Prior on length; a larger value results in longer')
            disp('                    tracts. If not specified, a suggested value will be')
            disp('                    computed based on the input data.')
            disp('order:              Number of polynomial coefficients to be optimized for')
            disp('                    each angle representing the tangent to the curve')
            disp('                    (default: 3).')
            disp('resH:               Discretization resolution (step size) of the Hough')
            disp('                    transform parameter space (default: 0.4).')
            disp('resS:               Discretization resolution of the tract arc length in')
            disp('                    the [-1,1]^3 box (default: 0.01).')
            disp('maxS:               Maximum partial tract length (in each side of the')
            disp('                    seed point) in the [-1,1]^3 box (default: 1).')
            disp('indexGPU:           If provided, the GPU device specified by this index')
            disp('                    will be used (default: 0 for no GPU).')
            disp('                    Run ''CSAODF_CLI GPUs'' to see the indices of the')
            disp('                    available GPUs. If N>1 GPUs are available, the user')
            disp('                    can specify a 1xN vector of GPU indices (e.g.,')
            disp('                    [1,2,3]) for ''indexGPU'', in addition to N for')
            disp('                    ''nCores'', so each CPU core uses a different GPU.')
            disp('nCores:             Number of CPU cores to be used. Choose Inf to use all')
            disp('                    the available cores (default: Inf if no GPU is used;')
            disp('                    1 if a GPU is used).')
            disp('createJobs:         If provided, ''numJobs'' batch job commands are')
            disp('                    generated and saved into ''jobFilename'' (which needs')
            disp('                    to be mandatorily specified). The jobs search')
            disp('                    different regions of the Hough-transform parameter')
            disp('                    space, and can be distributed and run in parallel,')
            disp('                    after which their results can be merged using')
            disp('                    ''CSAODF_CLI collectHoughResults''.')
            disp('jobFilename:        Name of the text file containing the batch job')
            disp('                    commands. See the above description for the')
            disp('                    ''createJobs'' option.')
            disp('verbose:            If more than 0, an estimate of the remaining time')
            disp('                    will be displayed every ''numIters'' iterations')
            disp('                    (default: 1000).')
            disp('nParts:             If provided, the Hough transform parameter space will')
            disp('                    be divided into ''nParts'' parts, and the code will')
            disp('                    only search the part number indicated by ''nPart'',')
            disp('                    which needs to be mandatorily specified (default: 1).')
            disp('                    Consider using the ''createJobs'' option instead.')
            disp('nPart:              See the description for ''nParts'' above (default: 1).')
            disp(' ')
            disp('Reference:')
            disp('I. Aganj, C. Lenglet, N. Jahanshad, E. Yacoub, N. Harel, P. Thompson, and')
            disp('G. Sapiro, "A Hough transform global probabilistic approach to')
            disp('multiple-subject diffusion MRI tractography," Medical Image Analysis,')
            disp('vol. 15, no. 4, pp. 414-425, 2011.')
            return
        end
        p.addRequired('input_ODFs');
        p.addRequired('input_seedPoints');
        p.addRequired('output_tracts');
        p.addRequired('mask');
        p.addParameter('pVol', []);
        p.addParameter('lambda', []);
        p.addParameter('nCores', []);
        p.addParameter('indexGPU', '0');
        p.addParameter('order', []);
        p.addParameter('resH', []);
        p.addParameter('resS', []);
        p.addParameter('maxS', []);
        p.addParameter('nParts', []);
        p.addParameter('nPart', []);
        p.addParameter('createJobs', []);
        p.addParameter('jobFilename', []);
        p.addParameter('verbose', []);
        p.parse(varargin{2:end});
        if isempty(p.Results.createJobs)
            if ~isempty(p.Results.jobFilename)
                error('jobFilename can only be provided together with createJobs.')
            end
            sh = load_nifti_allOS(p.Results.input_ODFs);
            mask = load_nifti_allOS(p.Results.mask);
            if isempty(p.Results.pVol)
                pVol = makeGFA(sh.vol);
            else
                pVol = load_nifti_allOS(p.Results.pVol);
                pVol = pVol.vol;
            end
            pVol = pVol.*(mask.vol>0);
            clear mask
            load(p.Results.input_seedPoints)
            param.indexGPU = str2num(p.Results.indexGPU);
            if isempty(p.Results.nCores)
                if ~any(param.indexGPU)
                    param.nCores = inf;
                else
                    param.nCores = length(param.indexGPU);
                end
            else
                param.nCores = str2double(p.Results.nCores);
            end
            if ~isempty(p.Results.lambda)
                param.lambda = str2double(p.Results.lambda);
            end
            if ~isempty(p.Results.order)
                param.order = str2double(p.Results.order);
            end
            if ~isempty(p.Results.resH)
                param.resH = str2double(p.Results.resH);
            end
            if ~isempty(p.Results.resS)
                param.resS = str2double(p.Results.resS);
            end
            if ~isempty(p.Results.maxS)
                param.maxS = str2double(p.Results.maxS);
            end
            if ~isempty(p.Results.nParts)
                param.nParts = str2double(p.Results.nParts);
            end
            if ~isempty(p.Results.nPart)
                param.nPart = str2double(p.Results.nPart);
            end
            if ~isempty(p.Results.verbose)
                param.verbose = str2double(p.Results.verbose);
            end
            tracts = HoughTract(sh.vol, pVol, seedPoints, box2vox, param);
            save(p.Results.output_tracts, 'tracts')
        else
            if ~isempty(p.Results.nParts) || ~isempty(p.Results.nPart)
                error('When creating parallel jobs, ''nParts'' and ''nPart'' should not be specified.')
            end
            if isempty(p.Results.jobFilename)
                error('When creating parallel jobs, jobFilename needs to be specified.')
            end
            nParts = str2double(p.Results.createJobs);
            fid = fopen(p.Results.jobFilename, 'w');
            for nPart = 1:nParts
                jobStr = [mfilename ' HoughTract ' p.Results.input_ODFs ' ' p.Results.input_seedPoints ' Part' num2str(nPart) 'of' num2str(nParts) '_' p.Results.output_tracts ' ' p.Results.mask ' nParts ' num2str(nParts) ' nPart ' num2str(nPart) ' indexGPU ' p.Results.indexGPU];
                if ~isempty(p.Results.nCores)
                    jobStr = [jobStr ' nCores ' p.Results.nCores];
                end
                if ~isempty(p.Results.pVol)
                    jobStr = [jobStr ' pVol ' p.Results.pVol];
                end
                if ~isempty(p.Results.lambda)
                    jobStr = [jobStr ' lambda ' p.Results.lambda];
                end
                if ~isempty(p.Results.order)
                    jobStr = [jobStr ' order ' p.Results.order];
                end
                if ~isempty(p.Results.resH)
                    jobStr = [jobStr ' resH ' p.Results.resH];
                end
                if ~isempty(p.Results.resS)
                    jobStr = [jobStr ' resS ' p.Results.resS];
                end
                if ~isempty(p.Results.maxS)
                    jobStr = [jobStr ' maxS ' p.Results.maxS];
                end
                if ~isempty(p.Results.verbose)
                    jobStr = [jobStr ' verbose ' p.Results.verbose];
                end
                fprintf(fid, [jobStr '\n']);
            end
            fclose(fid);
            disp([p.Results.jobFilename ' was created. To collect the job results, run:'])
            disp([mfilename ' collectHoughResults ' p.Results.output_tracts ' ' num2str(nParts)])
        end
    case 'COLLECTHOUGHRESULTS'
        if nargin<3
            disp('Collects and merges partial tractography results of batch jobs created by')
            disp('the ''createJobs'' option of ''CSAODF_CLI HoughTract''.')
            disp(' ')
            disp('CSAODF_CLI collectHoughResults <output_tracts> <nParts>')
            disp(' ')
            disp('output_tracts:   Name of the .mat file specified for ''output_tracts'' when')
            disp('                 creating the jobs.')
            disp('nParts:          The value specified for ''nParts'' when creating the jobs.')
            return
        end
        p.addRequired('output_tracts');
        p.addRequired('nParts');
        p.parse(varargin{2:end});
        nParts = str2double(p.Results.nParts);
        nResults = 0;
        for nPart = 1:nParts
            filename = ['Part' num2str(nPart) 'of' num2str(nParts) '_' p.Results.output_tracts];
            if exist(filename, 'file')
                load(filename)
                nResults = nResults + 1;
                partialTracts(nResults) = tracts;
            end
        end
        if nResults==0
            error(['No job results were found. They should look like:  Part*of' num2str(nParts) '_' p.Results.output_tracts])
        end
        sz2 = size(partialTracts(1).polyCoeffs);
        polyCoeffsAll = zeros([sz2 nResults]);
        partialLengthsAll = zeros(sz2(1),2,nResults);
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
        tracts.nParts = 1;
        tracts.nPart = 1;
        save(p.Results.output_tracts, 'tracts')
        disp([num2str(nResults) ' files were combined and saved in ' p.Results.output_tracts '.'])
    case 'SHOWTRACTS'
        if nargin<2
            disp('Visualizes the computed tracts.')
            disp(' ')
            disp('CSAODF_CLI showTracts <input_tracts> [overlay <overlay>] [nTracts <nTracts>]')
            disp('           [mask <mask>] [threshScore <threshScore>] [showSeeds <showSeeds>]')
            disp(' ')
            disp('input_tracts:   Name of the .mat file containing the tractography results.')
            disp(' ')
            disp('Optional parameters:')
            disp('overlay:        Image (in NIFTI), over which the tracts will be')
            disp('                visualized (default: no overlay).')
            disp('nTracts:        Number of top-scored tracts to be shown.')
            disp('mask:           Mask (in NIFTI) containing a region of interest. If')
            disp('                provided, only the tracts that intersect the mask will be')
            disp('                shown.')
            disp('threshScore:    If provided, only the tracts with a score above')
            disp('                ''threshScore'' will be shown.')
            disp('showSeeds:      If 1, seed points will be drawn as circles (default: 0).')
            return
        end
        p.addRequired('input_tracts');
        p.addParameter('overlay', []);
        p.addParameter('mask', []);
        p.addParameter('nTracts', []);
        p.addParameter('threshScore', []);
        p.addParameter('showSeeds', '0');
        p.parse(varargin{2:end});
        load(p.Results.input_tracts)
        if isempty(p.Results.overlay)
            overlayVol = [];
        else
            overlayVol = load_nifti_allOS(p.Results.overlay);
            overlayVol = overlayVol.vol;
        end
        if isempty(p.Results.mask)
            mask = [];
        else
            mask = load_nifti_allOS(p.Results.mask);
            mask = mask.vol;
        end
        if isempty(p.Results.nTracts)
            nTracts = [];
        else
            nTracts = str2double(p.Results.nTracts);
        end
        if isempty(p.Results.threshScore)
            threshScore = [];
        else
            threshScore = str2double(p.Results.threshScore);
        end
        showTracts(tracts, overlayVol, nTracts, mask, threshScore, str2double(p.Results.showSeeds))
    case 'EXPORT'
        if nargin<4
            disp('Exports the computed tracts to a .trk or text file.')
            disp('- The .trk file can be visualized by TrackVis (www.trackvis.org) and')
            disp('Freeview (https://surfer.nmr.mgh.harvard.edu).')
            disp('- The text file will include first the number of tracts, and then for each')
            disp('tract the number of points, the tract score, and the list of voxel')
            disp('coordinates of the points.')
            disp(' ')
            disp('CSAODF_CLI export <export_type> <input_tracts> <output_file>')
            disp('           [nTracts <nTracts>] [mask <mask>]')
            disp(' ')
            disp('export_type:    Can be either "trk" or "text".')
            disp('input_tracts:   The .mat file containing the tractography results.')
            disp('output_file:    Name of the output .trk or text file, where tracts')
            disp('                (sorted, high-score first) are exported.')
            disp(' ')
            disp('Optional parameters:')
            disp('nTracts:        Number of top-scored tracts to be exported.')
            disp('mask:           Mask (in NIFTI) containing a region of interest. If')
            disp('                provided, only the tracts that intersect the mask will be')
            disp('                exported.')
            return
        end
        p.addRequired('export_type', @(x) any(validatestring(x,{'trk', 'text'})));
        p.addRequired('input_tracts');
        p.addRequired('output_file');
        p.addParameter('mask', []);
        p.addParameter('nTracts', []);
        p.parse(varargin{2:end});
        load(p.Results.input_tracts)
        if isempty(p.Results.mask)
            mask = [];
        else
            mask = load_nifti_allOS(p.Results.mask);
            mask = mask.vol;
        end
        if isempty(p.Results.nTracts)
            nTracts = [];
        else
            nTracts = str2double(p.Results.nTracts);
        end
        switch validatestring(p.Results.export_type,{'trk', 'text'})
            case 'trk'
                exportTrackVis(p.Results.output_file, tracts, nTracts, mask);
            case 'text'
                exportText(p.Results.output_file, tracts, nTracts, mask);
        end
    case 'GPUS'
        nGPUs = gpuDeviceCount;
        switch nGPUs
            case 0
                disp('No GPU device was found. If this is unexpected, make sure that you have the latest CUDA driver.');
            case 1
                disp('The following GPU device was found.');
            otherwise
                disp(['The following ' num2str(nGPUs) ' GPU devices were found.'])
        end
        for n = 1:nGPUs
            fprintf('GPU device %d:\n\n', n)
            disp(gpuDevice(n))
        end
    otherwise
        disp('Usage:')
        disp('                      CSAODF_CLI <command> [options]')
        disp(' ')
        disp('One of the following can be used for ''command'':')
        disp(' ')
        disp('reconCSAODF:          Computes the q-ball CSA-ODFs.')
        disp('showODFs:             Visualizes the ODFs.')
        disp('verify:               Verifies the alignment of the gradient table with the image.')
        disp('makeGFA:              Computes the generalized fractional anisotropy.')
        disp('genSeedPoints:        Generates quasi-random seed points for tractography.')
        disp('HoughTract:           Performs the Hough-transform based tractography.')
        disp('collectHoughResults:  Merges partial tractography results of batch jobs.')
        disp('showTracts:           Visualizes the computed tracts.')
        disp('export:               Exports the computed tracts to a .trk or text file.')
        disp('GPUs:                 Shows the indices of the available GPUs.')
        disp(' ')
        disp('Developed by Iman Aganj.')
        disp('http://iman.mgh.harvard.edu')
        disp(' ')
end

end


function V = load_nifti_allOS(filename)

if ~isunix && strcmp(filename(end-2:end), '.gz')
    disp(['Extracting ' filename ' ...'])
    gunzip(filename)
    filename = filename(1:end-3);
end
V = load_nifti(filename);  % load_nifti.m is part of FreeSurfer's Matlab codes.

end


function save_nifti_allOS(V, filename)

if ~isunix && strcmp(filename(end-2:end), '.gz')
    filename = filename(1:end-3);
    save_nifti(V, filename);  % save_nifti.m is part of FreeSurfer's Matlab codes.
    disp(['Compressing ' filename ' ...'])
    gzip(filename)
else
    save_nifti(V, filename);
end

end
