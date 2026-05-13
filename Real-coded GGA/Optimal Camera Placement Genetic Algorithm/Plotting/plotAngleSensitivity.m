function plotAngleSensitivity(varargin)
% PLOTANGLESENSITIVITY  Sensitivity of total cost to angular perturbations
% of a single camera in the GA-best configuration.
%
% =====================================================================
% EXAMINER REVIEW
% =====================================================================
% What this plot claims to show
%   "Once the GA finds a configuration, how sharply does the cost rise
%   if a single camera is rotated about its three Euler axes (roll,
%   pitch, yaw)? Three subplots, one per axis, sweep the chosen
%   camera's angle from the optimum and re-evaluate the combined cost
%   function at each step. A narrow basin = the camera's orientation
%   is tightly tuned; a wide basin = orientation slack."
%
% Why this is useful for the thesis
%   - Speaks to physical realisability: if cost is flat over ±10°,
%     the lab technician does not need precision mounts.
%   - Identifies the angle the optimiser is most sensitive to —
%     useful argument when defending the cost-function design.
%   - One camera per call keeps each figure interpretable; loop
%     externally to cover all cameras if desired.
%
% USAGE
%   plotAngleSensitivity('TargetType', 1, 'CameraIndex', 3)
%   plotAngleSensitivity('TargetType', 2, 'CameraIndex', 1, ...
%                        'SweepRange', 60, 'SweepStep', 2)
%
% Name-Value Parameters
%   'TargetType', 'GridMode', 'Spacing', 'CostFunction', 'NumCameras',
%   'LogFile'      - same semantics as plotHeatmap_GAvsOptiTrack.
%   'CameraIndex'  - Which camera to perturb (1..NumCameras). Default 1.
%   'SweepRange'   - Half-range of the sweep in degrees. Default 45.
%   'SweepStep'    - Step size in degrees. Default 3.
%   'ShowAllCams'  - true to overlay faint sweep curves for every camera
%                    in the rig with the chosen camera highlighted.
%                    Default false.
%   'SaveAs'       - Output filename without extension.

    defaultLog = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                          'Results', 'Logs', 'GGA_RunsLog.mat');

    p = inputParser;
    addParameter(p, 'TargetType',   1,           @isnumeric);
    addParameter(p, 'GridMode',     1,           @isnumeric);
    addParameter(p, 'Spacing',      1.0,         @isnumeric);
    addParameter(p, 'CostFunction', 3,           @isnumeric);
    addParameter(p, 'NumCameras',   7,           @isnumeric);
    addParameter(p, 'LogFile',      defaultLog,  @ischar);
    addParameter(p, 'CameraIndex',  1,           @isnumeric);
    addParameter(p, 'SweepRange',   45,          @isnumeric);
    addParameter(p, 'SweepStep',    3,           @isnumeric);
    addParameter(p, 'ShowAllCams',  false,       @islogical);
    addParameter(p, 'SaveAs',       '',          @ischar);
    parse(p, varargin{:});
    opts = p.Results;

    sty = gaPlotStyle();
    ttStr = sty.TargetNames{opts.TargetType};
    gmStr = sty.GridNames{opts.GridMode};

    %% Load best GA run
    [gaChrom, specs, gaCost] = loadBestGARun(opts);

    numCams = specs.Cams;
    if opts.CameraIndex < 1 || opts.CameraIndex > numCams
        error('CameraIndex %d out of range 1..%d.', opts.CameraIndex, numCams);
    end

    %% Build the sweep grid (degrees -> radians)
    deltas_deg = -opts.SweepRange:opts.SweepStep:opts.SweepRange;
    nSteps     = numel(deltas_deg);

    %% Which cost-fn handle to evaluate
    cfHandle = pickCostHandle(opts.CostFunction);

    %% Pick cameras to sweep
    if opts.ShowAllCams
        camsToSweep = 1:numCams;
    else
        camsToSweep = opts.CameraIndex;
    end

    %% Storage: costs(c, axis, step), axis = 1 roll, 2 pitch, 3 yaw
    nSweep = numel(camsToSweep);
    costs  = nan(nSweep, 3, nSteps);

    fprintf(['plotAngleSensitivity: sweeping %d camera(s), %d steps per ' ...
             'axis (range ±%d°, step %d°).\n'], nSweep, nSteps, ...
             opts.SweepRange, opts.SweepStep);

    for cIdx = 1:nSweep
        cam = camsToSweep(cIdx);
        baseChrom = gaChrom;
        baseAng   = baseChrom((cam-1)*6 + (4:6));   % [alpha beta gamma] in rad

        for axis = 1:3
            for s = 1:nSteps
                dAng = deg2rad(deltas_deg(s));
                trial = baseChrom;
                trial((cam-1)*6 + 3 + axis) = baseAng(axis) + dAng;
                costs(cIdx, axis, s) = cfHandle(trial, specs);
            end
        end
        fprintf('  Camera %d done.\n', cam);
    end

    %% Plot — 1x3 subplots
    fig = figure('Name', sprintf('Angle sensitivity — %s, cam %d', ...
                                 ttStr, opts.CameraIndex), ...
        'Units', 'inches', ...
        'Position', [0.5, 0.5, sty.FigWidthDouble, sty.FigHeight + 0.4], ...
        'PaperPositionMode', 'auto', ...
        'Color', sty.BackgroundColor);

    axisNames  = {'Roll (\alpha)', 'Pitch (\beta)', 'Yaw (\gamma)'};
    chosenCol  = sty.CostFuncColors(3,:);    % green for the highlighted camera
    otherCol   = [0.55 0.55 0.55];           % grey for other cameras

    for axis = 1:3
        axSub = subplot(1, 3, axis);
        hold(axSub, 'on');

        legendH = gobjects(0);
        legendL = {};

        if opts.ShowAllCams
            for cIdx = 1:nSweep
                cam = camsToSweep(cIdx);
                y = squeeze(costs(cIdx, axis, :));
                if cam == opts.CameraIndex
                    h = plot(axSub, deltas_deg, y, '-', ...
                        'Color', chosenCol, 'LineWidth', sty.LineWidth + 0.4);
                    legendH(end+1) = h;                              %#ok<AGROW>
                    legendL{end+1} = sprintf('cam %d (chosen)', cam);
                else
                    plot(axSub, deltas_deg, y, '-', ...
                        'Color', [otherCol 0.45], ...
                        'LineWidth', sty.LineWidthThin, ...
                        'HandleVisibility', 'off');
                end
            end
            % Add a single grey handle for "other cameras" legend entry
            hG = plot(axSub, NaN, NaN, '-', ...
                'Color', otherCol, 'LineWidth', sty.LineWidthThin);
            legendH(end+1) = hG;                                     %#ok<AGROW>
            legendL{end+1} = 'other cams';
        else
            y = squeeze(costs(1, axis, :));
            h = plot(axSub, deltas_deg, y, '-', ...
                'Color', chosenCol, 'LineWidth', sty.LineWidth + 0.4);
            legendH(end+1) = h;                                      %#ok<AGROW>
            legendL{end+1} = sprintf('cam %d', opts.CameraIndex);
        end

        % Marker + reference line at optimum (delta = 0)
        plot(axSub, 0, gaCost, 'o', ...
            'MarkerFaceColor', chosenCol, 'MarkerEdgeColor', 'k', ...
            'MarkerSize', sty.MarkerSize, 'HandleVisibility', 'off');
        yL = ylim(axSub);
        plot(axSub, [0 0], yL, '--', ...
            'Color', [0.30 0.30 0.30 0.6], 'LineWidth', 0.8, ...
            'HandleVisibility', 'off');

        xlabel(axSub, sprintf('\\Delta%s (°)', axisNames{axis}), ...
            'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
        ylabel(axSub, 'Combined cost (CF3)', ...
            'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
        title(axSub, axisNames{axis}, ...
            'FontWeight', 'normal', ...
            'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
        set(axSub, 'FontSize', sty.FontSizeTick, 'FontName', sty.FontName, ...
            'Box', 'on', 'TickDir', 'out');
        grid(axSub, 'on');
        xlim(axSub, [-opts.SweepRange opts.SweepRange]);

        if axis == 3
            legend(axSub, legendH, legendL, ...
                'Location', 'best', 'FontSize', sty.FontSizeLegend);
        end
        hold(axSub, 'off');
    end

    sgtitle(sprintf('%s — Angle sensitivity of camera %d (GA-best, %s, sp=%.2f m, %d cams)', ...
        ttStr, opts.CameraIndex, gmStr, opts.Spacing, opts.NumCameras), ...
        'FontSize', sty.FontSizeTitle, 'FontWeight', 'bold', ...
        'FontName', sty.FontName);

    applyThesisStyle(fig);

    %% Export
    if isempty(opts.SaveAs)
        outName = sprintf('AngleSensitivity_%s_cam%d_%dC_GM%d_sp%.0fcm', ...
            ttStr, opts.CameraIndex, opts.NumCameras, opts.GridMode, ...
            opts.Spacing*100);
    else
        outName = opts.SaveAs;
    end
    exportgraphics(fig, [outName '.pdf'], ...
        'ContentType',     'vector', ...
        'BackgroundColor', sty.ExportBgColor);
    fprintf('Saved: %s.pdf\n', outName);
end


%% ---- Local helpers ------------------------------------------------------

function fn = pickCostHandle(cfType)
    switch cfType
        case 1
            fn = @resUncertaintyCost;
        case 2
            fn = @dynamicOcclusionCost;
        case 3
            fn = @combinedCostFunction;
        otherwise
            error('Unknown CostFunction %d (expected 1, 2, or 3).', cfType);
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
