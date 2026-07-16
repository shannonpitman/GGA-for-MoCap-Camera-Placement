function normTable = buildNormalisationSchedule(varargin)
% Table for normalisation calibration per instance for the combined (CF3)
% cost function

    addProjectPaths();

    p = inputParser;
    addParameter(p, 'NumCameras', [6 7 8], @isnumeric);
    addParameter(p, 'TargetTypes', [1 2], @isnumeric);
    addParameter(p, 'GridModes', [1 2], @isnumeric);
    addParameter(p, 'Spacing', 1.0, @isnumeric);
    addParameter(p, 'Replicates', 5, @isnumeric);
    addParameter(p, 'QuickGenerations', 50, @isnumeric);
    addParameter(p, 'QuickPopSize', 400, @isnumeric);
    addParameter(p, 'CamLowerBounds', [-5 -4.5 0   -pi -pi/2 -pi], @isnumeric);
    addParameter(p, 'CamUpperBounds', [ 5  4.5 4.8  pi  pi/2  pi], @isnumeric);
    addParameter(p, 'OutFile', '', @ischar);
    addParameter(p, 'Append',  false, @islogical); %if true, merge into prev file
    parse(p, varargin{:});
    opts = p.Results;

    projectRoot = fileparts(mfilename('fullpath'));
    if isempty(opts.OutFile)
        opts.OutFile = fullfile(projectRoot, 'Results', 'normTable.mat');
    end
    if ~isfolder(fileparts(opts.OutFile))
        mkdir(fileparts(opts.OutFile));
    end

    ttNames = {'UAV', 'UGV'};
    gmNames = {'Uniform', 'Normal'};

    % UGV target-slab settings (mirrors runCameraOptimiser.m).
    UGV_maxHeight = 0.5;
    UGV_zSpacing  = 0.25;

    % Suppress any incidental figures during the sweep.
    prevVis = get(0, 'DefaultFigureVisible');
    set(0, 'DefaultFigureVisible', 'off');
    cleanupObj = onCleanup(@() set(0, 'DefaultFigureVisible', prevVis));

    % Optionally resume from an existing file.
    existing = struct([]);
    if opts.Append && isfile(opts.OutFile)
        S = load(opts.OutFile, 'normTable');
        if isfield(S, 'normTable')
            existing = S.normTable;
        end
    end

    % Enumerate the instance grid.
    [TT, GM, NC] = ndgrid(opts.TargetTypes, opts.GridModes, opts.NumCameras);
    instances = [TT(:), GM(:), NC(:)];
    nInst = size(instances, 1);

    fprintf('\n Normalisation schedule: %d instance(s), %d replicate(s) each, %d gen x %d pop \n', nInst, opts.Replicates, opts.QuickGenerations, opts.QuickPopSize);

    normTable = struct([]);
    tAll = tic;

    for k = 1:nInst
        tt      = instances(k, 1);
        gm      = instances(k, 2);
        numCams = instances(k, 3);

        fprintf('\n[%d/%d] %s / %s | %d cams | %.2fm spacing\n', k, nInst, ttNames{tt}, gmNames{gm}, numCams, opts.Spacing);

        % Resume schedule handler
        if ~isempty(existing) && instanceExists(existing, tt, gm, numCams, opts.Spacing)
            fprintf('  already in %s -> skipped\n', opts.OutFile);
            row = existing(matchIdx(existing, tt, gm, numCams, opts.Spacing));
            normTable = appendRow(normTable, row);
            continue;
        end

        tInst = tic;

        % Build specs for this instance
        volume = [-4 4; -4 4; 0 4];
        if tt == 2
            volume(3, :)  = [0, UGV_maxHeight];
            targetSpacing = [opts.Spacing, opts.Spacing, min(UGV_zSpacing, UGV_maxHeight)];
        else
            targetSpacing = opts.Spacing;
        end

        specs = setupHardwareSpecs(numCams);
        specs.WeightUncertainty = 0.5;   % not used by CF1/CF2, set for consistency
        specs.WeightOcclusion = 0.5;
        specs.TargetType = tt;
        specs.TargetMode = gm;
        specs.Target = generateTargetSpace(volume, gm, targetSpacing);
        specs.NumPoints = size(specs.Target, 1);
        specs.spacing = opts.Spacing;
        if tt == 2
            specs.spacingZ = UGV_zSpacing;
        end
        specs.SectionCentres = generateSectionCentres(numCams, volume);
        specs.UseNormTable = false;   % calibration must not read its own output
        specs = setupCostParams(specs);
        specs.warmStart = false;
        specs.warmChromosomes = [];

        quickParams = setupGAparams(opts.QuickGenerations, opts.QuickPopSize);

        problem1 = setupProblem(numCams, 1, opts.CamUpperBounds, opts.CamLowerBounds);
        problem2 = setupProblem(numCams, 2, opts.CamUpperBounds, opts.CamLowerBounds);

        % Component 1: resolution (CF1), best of N replications
        fprintf('CF1 (resolution) x%d \n', opts.Replicates);
        [chromUnc, utopiaUnc, repUnc] = bestOfN(problem1, quickParams, specs, opts.Replicates);

        % Component 2: occlusion (CF2), best of N replications
        fprintf('CF2 (occlusion)  x%d \n', opts.Replicates);
        [chromOcc, utopiaOcc, repOcc] = bestOfN(problem2, quickParams, specs, opts.Replicates);

        % Pareto: opposite cost at each optimum
        nadirOcc = evalComponent(chromUnc, specs, 'occ');   % occlusion at the resolution optimum
        nadirUnc = evalComponent(chromOcc, specs, 'unc');   % resolution at the occlusion optimum

        uncertNorm = max(nadirUnc - utopiaUnc, eps);
        occlNorm   = max(nadirOcc - utopiaOcc, eps);

        % Record
        row = struct( ...
            'TargetType',   tt, ...
            'GridMode',     gm, ...
            'NumCameras',   numCams, ...
            'Spacing',      opts.Spacing, ...
            'utopiaUnc',    utopiaUnc, ...
            'nadirUnc',     nadirUnc, ...
            'uncertNorm',   uncertNorm, ...
            'utopiaOcc',    utopiaOcc, ...
            'nadirOcc',     nadirOcc, ...
            'occlNorm',     occlNorm, ...
            'chromUnc',     chromUnc, ...
            'chromOcc',     chromOcc, ...
            'replicateUnc', repUnc, ...
            'replicateOcc', repOcc, ...
            'elapsedSec',   toc(tInst));

        normTable = appendRow(normTable, row);

        fprintf('utopiaUnc=%.4g nadirUnc=%.4g -> uncertNorm=%.4g\n', utopiaUnc, nadirUnc, uncertNorm);
        fprintf('utopiaOcc=%.4g nadirOcc=%.4g -> occlNorm  =%.4g\n', utopiaOcc, nadirOcc, occlNorm);
        fprintf('(current setupCostParams constants: uncertNorm=%.4g occlNorm=%.4g)\n', specs.PreComputed.uncertNorm, specs.PreComputed.occlNorm);
        fprintf('instance done in %.1fs\n', toc(tInst));

        %Incremental save
        saveSchedule(opts.OutFile, normTable, opts, ttNames, gmNames);
    end

    saveSchedule(opts.OutFile, normTable, opts, ttNames, gmNames);
    fprintf('\n=== Schedule complete: %d instance(s) in %.1f min. Saved to %s ===\n', ...
        numel(normTable), toc(tAll) / 60, opts.OutFile);
end

%% Helpers
function [bestChrom, bestRaw, repCosts] = bestOfN(problem, params, specs, N)
    repCosts  = zeros(N, 1);
    bestRaw   = inf;
    bestChrom = [];
    for r = 1:N
        out = RunGA(problem, params, specs);
        repCosts(r) = out.bestsol.Cost;
        fprintf('run %d/%d: %.6g\n', r, N, out.bestsol.Cost);
        if out.bestsol.Cost < bestRaw
            bestRaw   = out.bestsol.Cost;
            bestChrom = out.bestsol.Chromosome;
        end
    end
end

function val = evalComponent(chrom, specs, which)
    [cameras, camCenters] = setupCameras(chrom, specs.Cams, specs.Resolution, specs.Focal, specs.FocalWide, specs.PrincipalPoint, specs.PixelSize);
    switch which
        case 'unc'
            val = resUncertainty(specs, cameras, camCenters);
        case 'occ'
            val = dynamicOcclusion(specs, cameras, camCenters);
        otherwise
            error('evalComponent:bad', 'which must be ''unc'' or ''occ''.');
    end
end

function tf = instanceExists(tbl, tt, gm, nc, sp)
    tf = ~isempty(matchIdx(tbl, tt, gm, nc, sp));
end

function idx = matchIdx(tbl, tt, gm, nc, sp)
    idx = find([tbl.TargetType] == tt & [tbl.GridMode] == gm & [tbl.NumCameras] == nc & abs([tbl.Spacing] - sp) < 1e-9, 1);
end

function tbl = appendRow(tbl, row)
    if isempty(tbl)
        tbl = row;
    else
        tbl(end + 1) = orderfields(row, tbl);
    end
end

function saveSchedule(outFile, normTable, opts, ttNames, gmNames)
% Persist the struct array, a readable table view, and provenance metadata.
    meta = struct( ...
        'CreatedOn',  datestr(now, 'yyyy-mm-dd HH:MM:SS'), ...
        'Method', 'payoff-table (utopia = diagonal, nadir = off-diagonal)', ...
        'Replicates', opts.Replicates, ...
        'QuickGenerations', opts.QuickGenerations, ...
        'QuickPopSize', opts.QuickPopSize, ...
        'Spacing', opts.Spacing, ...
        'CamLowerBounds', opts.CamLowerBounds, ...
        'CamUpperBounds', opts.CamUpperBounds);

    % Flat table view (scalars only) for quick inspection.
    if ~isempty(normTable)
        tt = [normTable.TargetType]';
        gm = [normTable.GridMode]';
        normSummary = table( ...
            tt, gm, [normTable.NumCameras]', [normTable.Spacing]', ...
            [normTable.utopiaUnc]', [normTable.nadirUnc]', [normTable.uncertNorm]', ...
            [normTable.utopiaOcc]', [normTable.nadirOcc]', [normTable.occlNorm]', ...
            'VariableNames', {'TargetType','GridMode','NumCameras','Spacing', ...
            'utopiaUnc','nadirUnc','uncertNorm','utopiaOcc','nadirOcc','occlNorm'});
        % ttNames is a row cell; indexing a vector follows the indexed
        % vector's orientation, so ttNames(tt) is 1xn -> transpose to nx1.
        normSummary.TargetTypeName = ttNames(tt)';
        normSummary.GridModeName   = gmNames(gm)';
    else
        normSummary = table();
    end

    save(outFile, 'normTable', 'normSummary', 'meta', '-v7.3');
end
