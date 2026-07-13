function calibrateNormalisation(varargin)

% USAGE:
% calibrateNormalisation() Default: UAV / Uniform, 7 cams, 1m
% calibrateNormalisation('TargetType', 2, 'GridMode', 2)   % UGV / Normal
% calibrateNormalisation('NumNadirSamples', 300, ...
%                           'QuickGenerations', 40, 'QuickPopSize', 60)
%
% Name-Value Parameters:
%   'NumCameras' - default 7
%   'TargetType' - 1 (UAV) or 2 (UGV). Default 1
%   'GridMode' - 1 (Uniform) or 2 (Normal). Default 1
%   'Spacing' - grid spacing in m. Default 1.0
%   'NumNadirSamples' - number of guided-random samples for the nadir estimate. Default 200
%   'QuickGenerations' - generations for the quick CF1/CF2 GA. Default 30
%   'QuickPopSize' - population for the quick CF1/CF2 GA. Default 50
%   'CamLowerBounds', 'CamUpperBounds' - as in batchRun7Cameras.m defaults

    addProjectPaths();

    p = inputParser;
    addParameter(p, 'NumCameras', 7, @isnumeric);
    addParameter(p, 'TargetType', 1, @isnumeric);
    addParameter(p, 'GridMode', 1, @isnumeric);
    addParameter(p, 'Spacing', 1.0, @isnumeric);
    addParameter(p, 'NumNadirSamples', 500,@isnumeric);
    addParameter(p, 'QuickGenerations', 75, @isnumeric);
    addParameter(p, 'QuickPopSize', 50, @isnumeric);
    addParameter(p, 'CamLowerBounds', [-5 -4.5 0   -pi -pi -pi], @isnumeric);
    addParameter(p, 'CamUpperBounds', [ 5  4.5 4.8  pi  pi  pi], @isnumeric);
    parse(p, varargin{:});
    opts = p.Results;

    ttNames = {'UAV', 'UGV'};
    gmNames = {'Uniform', 'Normal'};

    fprintf('\n Normalisation calibration: %s / %s, %d cams, %.1fm spacing \n',ttNames{opts.TargetType}, gmNames{opts.GridMode}, opts.NumCameras, opts.Spacing);

    %% Specs
    numCams = opts.NumCameras;
    volume  = [-4 4; -4 4; 0 4];
    if opts.TargetType == 2
        volume(3,:) = [0 0.5];
    end

    specs = setupHardwareSpecs(numCams);
    specs.WeightUncertainty = 0.5;
    specs.WeightOcclusion   = 0.5;
    specs.TargetType        = opts.TargetType;
    specs.TargetMode        = opts.GridMode;
    specs.Target            = generateTargetSpace(volume, opts.GridMode, opts.Spacing);
    specs.NumPoints         = size(specs.Target, 1);
    specs.spacing           = opts.Spacing;
    specs.SectionCentres    = generateSectionCentres(numCams, volume);
    specs                   = setupCostParams(specs);

    currentUncNorm  = specs.PreComputed.uncertNorm;
    currentOcclNorm = specs.PreComputed.occlNorm;

    %% Nadir: guided-random samples, no optimisation 
    VarMin = repmat(opts.CamLowerBounds, 1, numCams);
    VarMax = repmat(opts.CamUpperBounds, 1, numCams);

    fprintf('\n[1/2] Sampling %d guided-random chromosomes for the nadir estimate...\n', opts.NumNadirSamples);
    rawUnc = zeros(opts.NumNadirSamples, 1);
    rawOcc = zeros(opts.NumNadirSamples, 1);
    tSample = tic;
    for i = 1:opts.NumNadirSamples
        chrom = initialPopulation(VarMin, VarMax, specs.SectionCentres, numCams);
        [cameras, camCenters] = setupCameras(chrom, numCams, specs.Resolution, specs.Focal, specs.FocalWide, specs.PrincipalPoint, specs.PixelSize);
        rawUnc(i) = resUncertainty(specs, cameras, camCenters);
        rawOcc(i) = dynamicOcclusion(specs, cameras, camCenters);
    end
    fprintf('  done in %.1fs\n', toc(tSample));

    nadirUnc = mean(rawUnc);
    nadirOcc = mean(rawOcc);
    fprintf('  raw resUncertainty  : mean=%.4f  std=%.4f  range=[%.4f, %.4f]\n', ...
        nadirUnc, std(rawUnc), min(rawUnc), max(rawUnc));
    fprintf('  raw dynamicOcclusion: mean=%.4f  std=%.4f  range=[%.4f, %.4f]\n', ...
        nadirOcc, std(rawOcc), min(rawOcc), max(rawOcc));

    %% Utopia: single-objective GA runs 
    prevVis = get(0, 'DefaultFigureVisible');
    set(0, 'DefaultFigureVisible', 'off');
    cleanupObj = onCleanup(@() set(0, 'DefaultFigureVisible', prevVis)); 

    specs.warmStart       = false;
    specs.warmChromosomes = [];
    quickParams = setupGAparams(opts.QuickGenerations, opts.QuickPopSize);

    fprintf('\n[2/2] Quick CF1-only GA (%d gens x pop %d) -> utopia_res \n', ...
        opts.QuickGenerations, opts.QuickPopSize);
    problem1 = setupProblem(numCams, 1, opts.CamUpperBounds, opts.CamLowerBounds);
    tGA = tic;
    out1 = RunGA(problem1, quickParams, specs);
    utopiaUnc = out1.bestsol.Cost;
    fprintf('  done in %.1fs, utopia_res = %.4f\n', toc(tGA), utopiaUnc);

    fprintf('\nQuick CF2-only GA (%d gens x pop %d) -> utopia_occ \n', ...
        opts.QuickGenerations, opts.QuickPopSize);
    problem2 = setupProblem(numCams, 2, opts.CamUpperBounds, opts.CamLowerBounds);
    tGA = tic;
    out2 = RunGA(problem2, quickParams, specs);
    utopiaOcc = out2.bestsol.Cost;
    fprintf('  done in %.1fs, utopia_occ = %.4f\n', toc(tGA), utopiaOcc);

    %% Recommend new normalisation constants
    rangeUnc = max(nadirUnc - utopiaUnc, eps);
    rangeOcc = max(nadirOcc - utopiaOcc, eps);

    fprintf('\n--- Summary ---\n');
    fprintf('%-14s %-14s %-16s %-10s\n', 'Term', 'Utopia(best)', 'Nadir(typical)', 'Range');
    fprintf('%-14s %-14.4f %-16.4f %-10.4f\n', 'Resolution', utopiaUnc, nadirUnc, rangeUnc);
    fprintf('%-14s %-14.4f %-16.4f %-10.4f\n', 'Occlusion',  utopiaOcc, nadirOcc, rangeOcc);

    fprintf('\nCurrent constants   : uncertNorm = %-8.4f occlNorm = %-8.4f (ratio %.4f)\n', ...
        currentUncNorm, currentOcclNorm, currentUncNorm / currentOcclNorm);
    fprintf('Recommended (range) : uncertNorm = %-8.4f occlNorm = %-8.4f (ratio %.4f)\n', ...
        rangeUnc, rangeOcc, rangeUnc / rangeOcc);

    pctChangeUnc = 100 * (rangeUnc - currentUncNorm) / currentUncNorm;
    pctChangeOcc = 100 * (rangeOcc - currentOcclNorm) / currentOcclNorm;
    fprintf('  -> implies uncertNorm is %+.0f%% off, occlNorm is %+.0f%% off current values.\n', ...
        pctChangeUnc, pctChangeOcc);

    fprintf(['\nNOTE: report only. setupCostParams.m has NOT been changed.\n' ...
        'To see the effect on an already-saved GA-best chromosome, take the raw\n' ...
        'rawUnc/rawOcc values reportCostBreakdown.m already recovers and recompute:\n' ...
        '  newJres = w_res * rawUnc / %.4f;\n' ...
        '  newJocc = w_occ * rawOcc / %.4f;\n'], rangeUnc, rangeOcc);
end
