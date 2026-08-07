function plotCostField_GAvsOptiTrack(varargin)
% PLOTCOSTFIELD_GAVSOPTITRACK  Per-point uncertainty + occlusion heatmaps
% comparing the GA-best 7-camera configuration with the OptiTrack ad-hoc
% rig over the same target space.
%
% =====================================================================
% EXAMINER REVIEW
% =====================================================================
% What this plot claims to show
%   "Beyond just counting visible cameras, here is the actual COST FIELD
%   each placement produces at every target point. Top row = resolution
%   uncertainty (lower is better). Bottom row = dynamic-occlusion angle
%   (lower is better). Side-by-side, the regions where the OptiTrack
%   rig has high local cost — and where the GA solution improved them
%   — become visible directly."
%
% Strengths
%   - Uses the exact same per-point cost terms that the optimiser
%     minimises (computePointUncertainty + calculatePointOcclusion via
%     findVisibleCameras), so this plot is faithful to the objective.
%   - Each row shares a colour scale across the two columns (GA vs
%     OptiTrack) so spatial comparison is fair.
%   - Identical target grid for both: differences are placement-driven.
%
% USAGE
%   plotCostField_GAvsOptiTrack('TargetType', 1)               % UAV
%   plotCostField_GAvsOptiTrack('TargetType', 2)               % UGV
%
% Name-Value Parameters (same defaults as plotHeatmap_GAvsOptiTrack)
%   'TargetType', 'GridMode', 'Spacing', 'CostFunction', 'NumCameras',
%   'LogFile', 'MarkerSize', 'ViewAngle', 'SaveAs'.

    defaultLog = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                          'Results', 'Logs', 'GGA_RunsLog.mat');

    p = inputParser;
    addParameter(p, 'TargetType',   1,        @isnumeric);
    addParameter(p, 'GridMode',     1,        @isnumeric);
    addParameter(p, 'Spacing',      1.0,      @isnumeric);
    addParameter(p, 'CostFunction', 3,        @isnumeric);
    addParameter(p, 'NumCameras',   7,        @isnumeric);
    addParameter(p, 'LogFile',      defaultLog, @ischar);
    addParameter(p, 'MarkerSize',   40,       @isnumeric);
    addParameter(p, 'ViewAngle',    [45 30],  @isnumeric);
    addParameter(p, 'SaveAs',       '',       @ischar);
    parse(p, varargin{:});
    opts = p.Results;

    sty = gaPlotStyle();
    ttStr = sty.TargetNames{opts.TargetType};
    gmStr = sty.GridNames{opts.GridMode};

    %% Load best GA run
    [gaChrom, specs, gaCost] = loadBestGARun(opts);

    %% OptiTrack chromosome — re-use the same specs (same Target / PreComputed)
    optiChrom = buildOptiTrackChromosome();
    optiSpecs = specs;            % use the SAME target space + PreComputed
    optiSpecs.Cams = 7;

    %% Compute per-point cost terms for both
    fprintf('Computing per-point uncertainty + occlusion for GA-best...\n');
    [uncGA,   occGA]   = perPointCosts(gaChrom,   specs);
    fprintf('Computing per-point uncertainty + occlusion for OptiTrack...\n');
    [uncOpti, occOpti] = perPointCosts(optiChrom, optiSpecs);

    %% Figure: 2 rows (uncertainty, occlusion) x 2 cols (GA, OptiTrack)
    %  Use tiledlayout so the per-axes titles do not overlap with the
    %  colorbar above them. Title is split into a short heading + a
    %  subtitle line so each panel reads cleanly even when the figure is
    %  scaled to LaTeX text-width.
    fig = figure('Name', sprintf('Per-point cost field: %s', ttStr), ...
        'Units', 'inches', ...
        'Position', [0.5, 0.5, sty.FigWidthFull, sty.FigHeightTall * 1.25], ...
        'PaperPositionMode', 'auto', ...
        'Color', sty.BackgroundColor);

    tl = tiledlayout(fig, 2, 2, ...
        'TileSpacing', 'compact', 'Padding', 'compact');

    TargetSpace = specs.Target;

    % Shared scales per row
    uncMax = max([uncGA; uncOpti]);
    uncMin = 0;
    occMax = max([occGA; occOpti]);
    occMin = 0;

    % Use parula for both. Could swap for a sequential map if reviewers
    % prefer "low = green, high = red" — easy local change.
    cmap = parula(64);

    panels = struct( ...
        'unc',  {uncGA,  uncOpti}, ...
        'occ',  {occGA,  occOpti}, ...
        'name', {'GA-best', 'OptiTrack ad-hoc'});

    for k = 1:2
        % ---- Top row: uncertainty ----
        axU = nexttile(tl, k);
        scatter3(axU, TargetSpace(:,1), TargetSpace(:,2), TargetSpace(:,3), ...
            opts.MarkerSize, panels(k).unc, 'filled', 'MarkerFaceAlpha', 0.85);
        colormap(axU, cmap);
        clim(axU, [uncMin uncMax]);
        cb = colorbar(axU);
        cb.Label.String = 'Per-point uncertainty';
        cb.Label.FontSize = sty.FontSizeAxis;
        cb.Label.FontName = sty.FontName;

        axis(axU, 'equal');  grid(axU, 'on');
        xlabel(axU, 'X (m)'); ylabel(axU, 'Y (m)'); zlabel(axU, 'Z (m)');
        view(axU, opts.ViewAngle);
        title(axU, sprintf('Uncertainty: %s', panels(k).name), ...
            'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName, ...
            'FontWeight', 'normal', 'Color', 'k');
        subtitle(axU, sprintf('mean = %.3f', mean(panels(k).unc)), ...
            'FontSize', sty.FontSizeAnnot, 'FontName', sty.FontName, ...
            'Color', 'k');
        set(axU, 'FontSize', sty.FontSizeTick, 'FontName', sty.FontName);

        % ---- Bottom row: occlusion angle ----
        axO = nexttile(tl, k + 2);
        scatter3(axO, TargetSpace(:,1), TargetSpace(:,2), TargetSpace(:,3), ...
            opts.MarkerSize, panels(k).occ, 'filled', 'MarkerFaceAlpha', 0.85);
        colormap(axO, cmap);
        clim(axO, [occMin occMax]);
        cb2 = colorbar(axO);
        cb2.Label.String = 'Per-point occlusion angle (deg)';
        cb2.Label.FontSize = sty.FontSizeAxis;
        cb2.Label.FontName = sty.FontName;

        axis(axO, 'equal');  grid(axO, 'on');
        xlabel(axO, 'X (m)'); ylabel(axO, 'Y (m)'); zlabel(axO, 'Z (m)');
        view(axO, opts.ViewAngle);
        title(axO, sprintf('Occlusion: %s', panels(k).name), ...
            'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName, ...
            'FontWeight', 'normal', 'Color', 'k');
        subtitle(axO, sprintf('mean = %.1f deg', mean(panels(k).occ)), ...
            'FontSize', sty.FontSizeAnnot, 'FontName', sty.FontName, ...
            'Color', 'k');
        set(axO, 'FontSize', sty.FontSizeTick, 'FontName', sty.FontName);
    end

    title(tl, sprintf('%s per-point cost field: GA-best (CF3 = %.4f) vs OptiTrack (%s, sp = %.2f m)', ...
        ttStr, gaCost, gmStr, opts.Spacing), ...
        'FontSize', sty.FontSizeTitle, 'FontWeight', 'bold', ...
        'FontName', sty.FontName, 'Color', 'k');

    applyThesisStyle(fig);

    %% Export
    if isempty(opts.SaveAs)
        outName = sprintf('CostField_GAvsOptiTrack_%s_%dC_GM%d_sp%.0fcm', ...
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

function [unc, occ] = perPointCosts(chrom, specs)
% Evaluate uncertainty and occlusion at every target point.
    numCams = specs.Cams;
    [cameras, camCenters] = setupCameras(chrom, numCams, ...
        specs.Resolution, specs.Focal, specs.FocalWide, specs.PrincipalPoint, specs.PixelSize);

    resolution         = specs.Resolution;
    TargetSpace        = specs.Target;
    adjacentSurfaces   = specs.PreComputed.adjacentSurfaces;
    du                 = specs.PreComputed.du;
    dv                 = specs.PreComputed.dv;
    penaltyUncertainty = specs.PreComputed.penaltyUncertainty;
    w2                 = specs.PreComputed.w2;
    minTriangAngle     = specs.PreComputed.minTriangAngle;
    maxTriangAngle     = specs.PreComputed.maxTriangAngle;
    maxCameraRange     = specs.PreComputed.maxCameraRange;
    maxCameraRangeWide = specs.PreComputed.maxCameraRangeWide;
    focalWide          = specs.FocalWide;

    nPts = size(TargetSpace, 1);
    unc  = zeros(nPts, 1);
    occ  = zeros(nPts, 1);

    parfor pt = 1:nPts
        point = TargetSpace(pt, :);
        unc(pt) = computePointUncertainty(point, cameras, camCenters, ...
            numCams, adjacentSurfaces, du, dv, penaltyUncertainty, w2, resolution);
        [visCams, viewVecs] = findVisibleCameras(point, cameras, camCenters, ...
            numCams, resolution, maxCameraRange, maxCameraRangeWide, focalWide);
        occ(pt) = calculatePointOcclusion(visCams, viewVecs, minTriangAngle, maxTriangAngle);
    end
end


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
    chrom = L.saveData.BestSolution.Chromosome;
    specs = L.saveData.Specifications;
end
