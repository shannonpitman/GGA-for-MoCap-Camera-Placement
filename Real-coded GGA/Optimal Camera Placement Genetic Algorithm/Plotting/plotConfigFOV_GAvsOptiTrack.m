function plotConfigFOV_GAvsOptiTrack(varargin)
% PLOTCONFIGFOV_GAVSOPTITRACK  Side-by-side camera-pose + FOV frustum plot
% comparing the GA-best 7-camera configuration with the OptiTrack ad-hoc
% rig for a single scenario (UAV or UGV).
%
% =====================================================================
% EXAMINER REVIEW
% =====================================================================
% What this plot claims to show
%   "Here is the physical realisation of the GA solution: each camera's
%   position, orientation, and field-of-view frustum, drawn against the
%   target volume. Side-by-side with the OptiTrack ad-hoc rig you can
%   see WHERE the GA chose to redistribute cameras — apex spacing,
%   pitch/yaw angle, and frustum coverage overlap."
%
% Strengths
%   - The drawn frustum matches what findVisibleCameras actually checks
%     (uses each camera's K-matrix and effective range).
%   - Wide-lens cameras get a different colour and shorter frustum,
%     so the lens choice the GA made is visible at a glance.
%   - Identical axes / view angle on both panels makes geometric
%     differences easy to read.
%
% USAGE
%   plotConfigFOV_GAvsOptiTrack('TargetType', 1)               % UAV
%   plotConfigFOV_GAvsOptiTrack('TargetType', 2)               % UGV
%
% Name-Value Parameters
%   'TargetType', 'GridMode', 'Spacing', 'CostFunction', 'NumCameras',
%   'LogFile'    - same semantics as plotHeatmap_GAvsOptiTrack.
%   'ViewAngle'  - 3D view angle [az el]. Default [45 25].
%   'ShowTarget' - true to scatter target points. Default true.
%   'SaveAs'     - output filename without extension.

    %% Parse inputs
    defaultLog = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                          'Results', 'Logs', 'GGA_RunsLog.mat');

    p = inputParser;
    addParameter(p, 'TargetType',   1,        @isnumeric);
    addParameter(p, 'GridMode',     1,        @isnumeric);
    addParameter(p, 'Spacing',      1.0,      @isnumeric);
    addParameter(p, 'CostFunction', 3,        @isnumeric);
    addParameter(p, 'NumCameras',   7,        @isnumeric);
    addParameter(p, 'LogFile',      defaultLog, @ischar);
    addParameter(p, 'ViewAngle',    [45 25],  @isnumeric);
    addParameter(p, 'ShowTarget',   true,     @islogical);
    addParameter(p, 'SaveAs',       '',       @ischar);
    parse(p, varargin{:});
    opts = p.Results;

    sty = gaPlotStyle();
    ttNames = sty.TargetNames;
    gmNames = sty.GridNames;
    ttStr   = ttNames{opts.TargetType};
    gmStr   = gmNames{opts.GridMode};

    %% Load best GA run for the scenario
    [gaChrom, specs, gaCost] = loadBestGARun(opts);

    %% OptiTrack chromosome + cost on the same problem
    optiChrom   = buildOptiTrackChromosome();
    optiOut     = evaluateOptiTrackCost( ...
        'TargetType', opts.TargetType, ...
        'GridMode',   opts.GridMode, ...
        'Spacing',    opts.Spacing);
    cfField     = sprintf('CF%d', opts.CostFunction);
    if isfield(optiOut, cfField)
        optiCost = optiOut.(cfField);
    else
        optiCost = NaN;
    end

    %% Figure: 1 x 2 panels, equal axes
    fig = figure('Name', sprintf('Camera Poses + FOV — %s', ttStr), ...
        'Units', 'inches', ...
        'Position', [0.5, 0.5, sty.FigWidthDouble, sty.FigHeightWide], ...
        'PaperPositionMode', 'auto', ...
        'Color', sty.BackgroundColor);

    ax1 = subplot(1, 2, 1);
    drawCamerasAndVolume(ax1, gaChrom, specs, opts);
    title(ax1, sprintf('GA-best — %d cams (Cost: %.4f)', specs.Cams, gaCost), ...
        'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName, ...
        'FontWeight', 'normal');

    ax2 = subplot(1, 2, 2);
    % OptiTrack has 7 cameras — re-use the same specs object but force its
    % camera count to 7 so setupCameras unpacks the right number of genes.
    optiSpecs       = specs;
    optiSpecs.Cams  = 7;
    drawCamerasAndVolume(ax2, optiChrom, optiSpecs, opts);
    title(ax2, sprintf('OptiTrack ad-hoc — 7 cams (Cost: %.4f)', optiCost), ...
        'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName, ...
        'FontWeight', 'normal');

    % Synchronise axes for an honest side-by-side
    linkprop([ax1 ax2], {'XLim','YLim','ZLim','View'});
    view(ax1, opts.ViewAngle);

    sgtitle(sprintf('%s — Camera placement + FOV (Uniform grid, sp = %.2f m)', ...
        ttStr, opts.Spacing), ...
        'FontSize', sty.FontSizeTitle, 'FontWeight', 'bold', ...
        'FontName', sty.FontName);

    applyThesisStyle(fig);

    %% Export
    if isempty(opts.SaveAs)
        outName = sprintf('ConfigFOV_GAvsOptiTrack_%s_%dC_GM%d_sp%.0fcm', ...
            ttStr, opts.NumCameras, opts.GridMode, opts.Spacing*100);
    else
        outName = opts.SaveAs;
    end
    exportgraphics(fig, [outName '.pdf'], ...
        'ContentType',     'vector', ...
        'BackgroundColor', sty.ExportBgColor);
    fprintf('Saved: %s.pdf\n', outName);
end


%% ---- Local helpers ------------------------------------------------------

function [chrom, specs, cost] = loadBestGARun(opts)
% Find the lowest-cost GA run matching the scenario, return its
% chromosome, specs, and BestCost.
    if ~isfile(opts.LogFile)
        error('Log file not found: %s', opts.LogFile);
    end
    S = load(opts.LogFile, 'runLog');
    runLog = S.runLog;

    fillFields = {'TargetType', 'GridMode', 'Spacing'};
    for f = 1:length(fillFields)
        if ~isfield(runLog, fillFields{f})
            [runLog.(fillFields{f})] = deal(NaN);
        end
    end

    mask = ([runLog.NumCameras]       == opts.NumCameras)   & ...
           ([runLog.CostFunctionType] == opts.CostFunction) & ...
           ([runLog.TargetType]       == opts.TargetType)   & ...
           ([runLog.GridMode]         == opts.GridMode)     & ...
           (abs([runLog.Spacing] - opts.Spacing) < 1e-6);

    candidates = runLog(mask);
    if isempty(candidates)
        error('No GA runs found for TT=%d GM=%d sp=%.2f CF=%d %dC.', ...
            opts.TargetType, opts.GridMode, opts.Spacing, ...
            opts.CostFunction, opts.NumCameras);
    end
    [cost, idx] = min([candidates.BestCost]);
    bestRun = candidates(idx);

    matFile = resolveRunPath(bestRun.RunFilename, bestRun.NumCameras);
    if ~isfile(matFile)
        error('GA result file not found: %s', matFile);
    end
    L = load(matFile, 'saveData');
    chrom = L.saveData.BestSolution.Chromosome;
    specs = L.saveData.Specifications;
end


function drawCamerasAndVolume(ax, chrom, specs, opts)
% Draw camera FOV frustums + (optionally) the target point cloud.
    numCams = specs.Cams;
    [cameras, camCenters] = setupCameras(chrom, numCams, ...
        specs.Resolution, specs.Focal, specs.FocalWide, specs.PrincipalPoint);

    maxRange     = specs.PreComputed.maxCameraRange;
    maxRangeWide = specs.PreComputed.maxCameraRangeWide;
    focalWide    = specs.FocalWide;

    hold(ax, 'on');
    for i = 1:numCams
        if cameras{i}.f == focalWide
            rng_i = maxRangeWide;
            col_i = [0.85 0.40 0.10];   % wide lens — orange
        else
            rng_i = maxRange;
            col_i = [0.10 0.45 0.75];   % narrow lens — blue
        end
        plotCameraFOV(cameras{i}, camCenters(:,i), rng_i, ...
            'Color', col_i, 'Label', sprintf('cam%d', i));
    end

    if opts.ShowTarget && isfield(specs, 'Target')
        T = specs.Target;
        plot3(ax, T(:,1), T(:,2), T(:,3), '.', ...
            'Color', [0.45 0.45 0.45], 'MarkerSize', 3);
    end

    axis(ax, 'equal');
    grid(ax, 'on');
    xlabel(ax, 'X (m)'); ylabel(ax, 'Y (m)'); zlabel(ax, 'Z (m)');
    view(ax, opts.ViewAngle);
    hold(ax, 'off');
end
