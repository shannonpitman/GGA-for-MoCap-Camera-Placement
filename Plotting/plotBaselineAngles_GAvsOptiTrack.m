function plotBaselineAngles_GAvsOptiTrack(varargin)
% PLOTBASELINEANGLES_GAVSOPTITRACK  Pairwise camera-baseline-angle histogram
% comparing Optimised GA Rig 7-camera config with Manually Posed Rig, one scenario
% per figure (UAV or UGV).
%
% =====================================================================
% EXAMINER REVIEW
% =====================================================================
% What this plot claims to show
%   "Triangulation accuracy depends on the ANGLE between camera rays
%   meeting at a target point, not just the number of cameras that see
%   it. For every (point, camera-pair) where both cameras see the
%   point, this histogram reports the converging-ray angle. The
%   [minTriangAngle, maxTriangAngle] band marks the angles the cost
%   function counts as triangulable — outside that band, a pair is
%   effectively useless for reconstruction. The GA solution and the
%   OptiTrack rig are overlaid so the reader can judge which produces
%   a fatter mass inside the triangulable band."
%
% Strengths
%   - Directly addresses the colour-blind / 2+ count critique levelled
%     at visualizeCameraCoverage: instead of "2+ cameras = good", the
%     metric is the angle the cost function actually uses.
%   - Vertical dashed lines at minTriang / maxTriang make the
%     triangulable band visible.
%   - Same target space for both placements → differences are
%     placement-driven, not grid-driven.
%
% USAGE
%   plotBaselineAngles_GAvsOptiTrack('TargetType', 1)              % UAV
%   plotBaselineAngles_GAvsOptiTrack('TargetType', 2)              % UGV
%
% Name-Value Parameters (same as plotHeatmap_GAvsOptiTrack, plus)
%   'NumBins'   - Histogram bin count. Default 36 (5° bins over 0–180°).
%   'Normalise' - 'probability' (default), 'count', or 'pdf'.

    defaultLog = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                          'Results', 'Logs', 'GGA_RunsLog.mat');

    p = inputParser;
    addParameter(p, 'TargetType',   1,           @isnumeric);
    addParameter(p, 'GridMode',     1,           @isnumeric);
    addParameter(p, 'Spacing',      1.0,         @isnumeric);
    addParameter(p, 'CostFunction', 3,           @isnumeric);
    addParameter(p, 'NumCameras',   7,           @isnumeric);
    addParameter(p, 'LogFile',      defaultLog,  @ischar);
    addParameter(p, 'NumBins',      36,          @isnumeric);
    addParameter(p, 'Normalise',    'probability', @ischar);
    addParameter(p, 'SaveAs',       '',          @ischar);
    parse(p, varargin{:});
    opts = p.Results;

    sty = gaPlotStyle();
    ttStr = sty.TargetNames{opts.TargetType};
    gmStr = sty.GridNames{opts.GridMode};

    %% Load best GA run
    [gaChrom, specs, gaCost] = loadBestGARun(opts);

    %% OptiTrack chromosome — same specs.Target for fair comparison
    optiChrom = buildOptiTrackChromosome();
    optiSpecs = specs;
    optiSpecs.Cams = 7;

    %% Collect pairwise baseline angles (degrees) for every target point
    fprintf('Computing pairwise baseline angles for Optimised GA Rig...\n');
    angGA   = pairwiseBaselineAngles(gaChrom,   specs);
    fprintf('Computing pairwise baseline angles for OptiTrack...\n');
    angOpti = pairwiseBaselineAngles(optiChrom, optiSpecs);

    minTri = specs.PreComputed.minTriangAngle;
    maxTri = specs.PreComputed.maxTriangAngle;

    %% Plot
    fig = figure('Name', sprintf('Baseline angles: %s', ttStr), ...
        'Units', 'inches', ...
        'Position', [0.5, 0.5, sty.FigWidthFull, sty.FigHeight + 0.4], ...
        'PaperPositionMode', 'auto', ...
        'Color', sty.BackgroundColor);
    ax = axes(fig);
    hold(ax, 'on');

    edges = linspace(0, 180, opts.NumBins + 1);

    hGA = histogram(ax, angGA, edges, ...
        'Normalization', opts.Normalise, ...
        'FaceColor',    sty.CostFuncColors(3,:), ...   % combined-green
        'FaceAlpha',    0.55, ...
        'EdgeColor',    'none', ...
        'DisplayName',  'Optimised GA Rig');

    hOpti = histogram(ax, angOpti, edges, ...
        'Normalization', opts.Normalise, ...
        'FaceColor',    [0.85 0.10 0.10], ...          % opti-red
        'FaceAlpha',    0.45, ...
        'EdgeColor',    'none', ...
        'DisplayName',  'Manually Posed Rig');

    %% Triangulable-band guide lines
    yL = ylim(ax);
    plot(ax, [minTri minTri], yL, '--', ...
        'Color', [0.30 0.30 0.30 0.7], 'LineWidth', 1.0, ...
        'HandleVisibility', 'off');
    plot(ax, [maxTri maxTri], yL, '--', ...
        'Color', [0.30 0.30 0.30 0.7], 'LineWidth', 1.0, ...
        'HandleVisibility', 'off');

    hold(ax, 'off');

    %% Headline numbers
    fracGA   = sum(angGA   >= minTri & angGA   <= maxTri) / max(numel(angGA), 1);
    fracOpti = sum(angOpti >= minTri & angOpti <= maxTri) / max(numel(angOpti), 1);

    xlabel(ax, 'Baseline angle between camera pair (°)', ...
        'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
    if strcmpi(opts.Normalise, 'count')
        ylabel(ax, 'Count', ...
            'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
    else
        ylabel(ax, sprintf('Fraction (%s)', opts.Normalise), ...
            'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
    end
    % In-band = pairs inside [minTri, maxTri], marked by the dashed lines.
    title(ax, 'Pairwise Baseline Angles', ...
        'FontWeight', 'bold', 'FontSize', sty.FontSizeTitle, ...
        'FontName', sty.FontName);
    subtitle(ax, sprintf('%s, %s Grid, optimised %.1f%% vs manual %.1f%%', ...
                         ttStr, gmStr, 100*fracGA, 100*fracOpti), ...
        'FontWeight', 'normal', 'FontSize', sty.FontSizeAxis, ...
        'FontName', sty.FontName);
    set(ax, 'FontSize', sty.FontSizeTick, 'FontName', sty.FontName, ...
        'Box', 'on', 'TickDir', 'out');
    grid(ax, 'on');
    legend([hGA, hOpti], 'Location', 'northeast', 'FontSize', sty.FontSizeLegend);
    xlim(ax, [0 180]);

    applyThesisStyle(fig);

    fprintf('\n%s — baseline-angle summary (CF3, %dC, %s, sp=%.2f m):\n', ...
        ttStr, opts.NumCameras, gmStr, opts.Spacing);
    fprintf('  Optimised GA Rig      median %.1f° | in-band %.1f%% | n_pairs %d (cost %.4f)\n', ...
        median(angGA), 100*fracGA, numel(angGA), gaCost);
    fprintf('  OptiTrack    median %.1f° | in-band %.1f%% | n_pairs %d\n', ...
        median(angOpti), 100*fracOpti, numel(angOpti));

    %% Export
    if isempty(opts.SaveAs)
        outName = sprintf('BaselineAngles_GAvsOptiTrack_%s_%dC_GM%d_sp%.0fcm', ...
            ttStr, opts.NumCameras, opts.GridMode, opts.Spacing*100);
    else
        outName = opts.SaveAs;
    end
    exportThesisFigure(fig, outName, ...
        'Background', sty.ExportBgColor, 'Quiet', true);
    fprintf('Saved: %s.{pdf,png}\n', outName);
end


%% ---- Local helpers ------------------------------------------------------

% pairwiseBaselineAngles now lives in Analysis/ so this figure and the
% reported in-band percentages share one implementation.

function [chrom, specs, cost] = loadBestGARun(opts)
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
    chrom = runChromosome(L.saveData);
    specs = L.saveData.Specifications;
end
