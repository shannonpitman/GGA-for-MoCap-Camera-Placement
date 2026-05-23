function plotConfigFOV_GAvsOptiTrack(varargin)
% Side-by-side camera-pose plot comparing the GA-best 7-camera
% configuration with the OptiTrack ad-hoc arrangement for a single
% scenario (UAV or UGV). Each camera is drawn as a small fixed-size
% pyramid pointing along its optical axis, coloured by lens type
% (narrow vs wide). The full FOV frustum is no longer drawn because it
% cluttered the figure and obscured the pose comparison.
%
% USAGE
%   plotConfigFOV_GAvsOptiTrack('TargetType', 1)               % UAV
%   plotConfigFOV_GAvsOptiTrack('TargetType', 2)               % UGV
%
% Name-Value Parameters
%   'TargetType', 'GridMode', 'Spacing', 'CostFunction', 'NumCameras',
%   'LogFile'      - same semantics as plotHeatmap_GAvsOptiTrack.
%   'ViewAngle'    - 3D view angle [az el]. Default [45 25].
%   'ShowTarget'   - true to scatter target points. Default true.
%   'PyramidScale' - length (m) of each pose pyramid along its optical
%                    axis. Default 0.4 m. Base width = PyramidScale/2.
%   'SaveAs'       - output filename without extension.

    %% Parse inputs
    defaultLog = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'Results', 'Logs', 'GGA_RunsLog.mat');

    p = inputParser;
    addParameter(p, 'TargetType',   1,        @isnumeric);
    addParameter(p, 'GridMode',     1,        @isnumeric);
    addParameter(p, 'Spacing',      1.0,      @isnumeric);
    addParameter(p, 'CostFunction', 3,        @isnumeric);
    addParameter(p, 'NumCameras',   7,        @isnumeric);
    addParameter(p, 'LogFile',      defaultLog, @ischar);
    addParameter(p, 'ViewAngle',    [45 25],  @isnumeric);
    addParameter(p, 'ShowTarget',   true,     @islogical);
    addParameter(p, 'PyramidScale', 0.4,      @(x) isnumeric(x) && isscalar(x) && x > 0);
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
    %  Use tiledlayout with loose padding so the per-axes title and the
    %  northeast lens legend stop competing for the same band of pixels.
    %  Figure height bumped slightly for the same reason.
    fig = figure('Name', sprintf('Camera poses: %s', ttStr), ...
        'Units', 'inches', ...
        'Position', [0.5, 0.5, sty.FigWidthDouble, sty.FigHeightWide + 0.6], ...
        'PaperPositionMode', 'auto', ...
        'Color', sty.BackgroundColor);

    tl = tiledlayout(fig, 1, 2, ...
        'TileSpacing', 'compact', 'Padding', 'loose');

    ax1 = nexttile(tl);
    drawCamerasAndVolume(ax1, gaChrom, specs, opts);
    % Per-axes title kept short (no cost number) so the lens legend
    % at northeast does not run into the heading. Cost moved to sgtitle.
    title(ax1, sprintf('GA-best (%d cams)', specs.Cams), ...
        'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName, ...
        'FontWeight', 'normal', 'Color', 'k');
    addLensLegend(ax1, sty);

    ax2 = nexttile(tl);
    % OptiTrack has 7 cameras — re-use the same specs object but force its
    % camera count to 7 so setupCameras unpacks the right number of genes.
    optiSpecs       = specs;
    optiSpecs.Cams  = 7;
    drawCamerasAndVolume(ax2, optiChrom, optiSpecs, opts);
    title(ax2, 'OptiTrack ad-hoc (7 cams)', ...
        'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName, ...
        'FontWeight', 'normal', 'Color', 'k');
    addLensLegend(ax2, sty);

    % Synchronise axes for an honest side-by-side
    linkprop([ax1 ax2], {'XLim','YLim','ZLim','View'});
    view(ax1, opts.ViewAngle);

    title(tl, sprintf('%s camera placement (Uniform grid, sp = %.2f m): GA-best cost %.4f vs OptiTrack %.4f', ...
        ttStr, opts.Spacing, gaCost, optiCost), ...
        'FontSize', sty.FontSizeTitle, 'FontWeight', 'bold', ...
        'FontName', sty.FontName, 'Color', 'k');

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
% Draw room-volume wireframe + mini pose pyramids for each camera + (optionally) the target point cloud.
%
% Each camera is rendered as a small fixed-size pyramid whose apex is at
% the camera centre and whose base lies at apex + opticalAxis *
% PyramidScale. Colour encodes lens type:
%   narrow lens -> blue   (0.10 0.45 0.75)
%   wide   lens -> orange (0.85 0.40 0.10)
% A small filled marker at the apex draws the eye to the lens centre.

    numCams = specs.Cams;
    [cameras, camCenters] = setupCameras(chrom, numCams, ...
        specs.Resolution, specs.Focal, specs.FocalWide, specs.PrincipalPoint, specs.PixelSize);

    focalWide    = specs.FocalWide;
    pyrLen       = opts.PyramidScale;
    pyrHalfBase  = 0.5 * pyrLen / 2;       % base width = pyrLen/2

    narrowCol = [0.10 0.45 0.75];
    wideCol   = [0.85 0.40 0.10];

    hold(ax, 'on');

    % --- Room-volume wireframe (gives the cameras spatial reference) ---
    if isfield(specs, 'Target') && ~isempty(specs.Target)
        T = specs.Target;
        bb = [min(T,[],1); max(T,[],1)];   % 2x3 [x;y;z] bounds
        drawBoxWire(ax, bb, [0.30 0.30 0.30], 0.6, 0.9);
    end

    % --- Optional target point cloud (drawn first so cameras sit on top) ---
    if opts.ShowTarget && isfield(specs, 'Target')
        T = specs.Target;
        plot3(ax, T(:,1), T(:,2), T(:,3), '.', ...
            'Color', [0.55 0.55 0.55 0.45], 'MarkerSize', 3, ...
            'HandleVisibility', 'off');
    end

    % --- Mini pose pyramids ---
    for i = 1:numCams
        if cameras{i}.f == focalWide
            col_i = wideCol;
        else
            col_i = narrowCol;
        end
        R    = cameras{i}.T.rotm;        % camera-to-world rotation
        axisDir = R(:, 3);               % optical axis (camera +Z) in world
        rightV  = R(:, 1);
        upV     = R(:, 2);
        drawPosePyramid(ax, camCenters(:,i), axisDir, rightV, upV, ...
            pyrLen, pyrHalfBase, col_i);
        % Small dot at the apex so the lens centre is obvious.
        plot3(ax, camCenters(1,i), camCenters(2,i), camCenters(3,i), 'o', ...
            'MarkerSize', 5, 'MarkerFaceColor', col_i, ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.6, ...
            'HandleVisibility', 'off');
        % Small label next to apex
        text(ax, camCenters(1,i), camCenters(2,i), camCenters(3,i), ...
            sprintf('  cam%d', i), 'Color', col_i, ...
            'FontSize', 8, 'FontWeight', 'bold', ...
            'HandleVisibility', 'off');
    end

    axis(ax, 'equal');
    grid(ax, 'on');
    xlabel(ax, 'X (m)'); ylabel(ax, 'Y (m)'); zlabel(ax, 'Z (m)');
    view(ax, opts.ViewAngle);
    hold(ax, 'off');
end


function drawPosePyramid(ax, apex, axisDir, rightV, upV, len, halfBase, col)
% Draw a small 4-sided pyramid with apex at `apex` pointing along
% `axisDir`, with base square sized by `halfBase`. The base sits at
% apex + axisDir * len. Faces are drawn translucent; edges drawn opaque
% so the wireframe pyramid reads cleanly in 3D.
    apex     = apex(:);
    axisDir  = axisDir(:) / norm(axisDir);
    rightV   = rightV(:)  / norm(rightV);
    upV      = upV(:)     / norm(upV);
    baseCen  = apex + axisDir * len;

    % 4 base corners (square in the right/up plane at base centre)
    corners = [ baseCen + halfBase*rightV + halfBase*upV, ...
                baseCen - halfBase*rightV + halfBase*upV, ...
                baseCen - halfBase*rightV - halfBase*upV, ...
                baseCen + halfBase*rightV - halfBase*upV ];   % 3x4

    % Side faces (translucent)
    for c = 1:4
        c2 = mod(c, 4) + 1;
        triX = [apex(1), corners(1,c), corners(1,c2)];
        triY = [apex(2), corners(2,c), corners(2,c2)];
        triZ = [apex(3), corners(3,c), corners(3,c2)];
        patch(ax, triX, triY, triZ, col, ...
            'FaceAlpha', 0.30, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
    end

    % Base quad (slightly stronger alpha so the orientation reads)
    patch(ax, corners(1,:), corners(2,:), corners(3,:), col, ...
        'FaceAlpha', 0.45, 'EdgeColor', col, 'LineWidth', 1.0, ...
        'HandleVisibility', 'off');

    % Apex-to-corner edges
    for c = 1:4
        plot3(ax, [apex(1) corners(1,c)], ...
                  [apex(2) corners(2,c)], ...
                  [apex(3) corners(3,c)], ...
              '-', 'Color', col, 'LineWidth', 1.0, ...
              'HandleVisibility', 'off');
    end
end


function addLensLegend(ax, sty)
% Add a per-subplot legend in the top-right of the axes with two
% entries: Narrow lens (blue) and Wide lens (orange). The legend
% entries are invisible NaN-data scatter plots so they appear in the
% legend without polluting the 3D scene.
    narrowCol = [0.10 0.45 0.75];
    wideCol   = [0.85 0.40 0.10];
    hold(ax, 'on');
    hN = plot3(ax, NaN, NaN, NaN, 's', ...
        'MarkerFaceColor', narrowCol, 'MarkerEdgeColor', 'k', ...
        'MarkerSize', 9, 'LineStyle', 'none', ...
        'DisplayName', 'Narrow lens');
    hW = plot3(ax, NaN, NaN, NaN, 's', ...
        'MarkerFaceColor', wideCol,   'MarkerEdgeColor', 'k', ...
        'MarkerSize', 9, 'LineStyle', 'none', ...
        'DisplayName', 'Wide lens');
    lgd = legend(ax, [hN hW], {'Narrow lens', 'Wide lens'}, ...
        'Location', 'northeast', ...
        'FontSize', sty.FontSizeLegend);
    set(lgd, 'Color', 'w', 'EdgeColor', 'k', 'Box', 'on', 'TextColor', 'k');
    hold(ax, 'off');
end


function drawBoxWire(ax, bb, col, lw, alpha)
% Draw the 12 edges of an axis-aligned bounding box bb = [xmin ymin zmin;
% xmax ymax zmax].
    x = bb(:,1);  y = bb(:,2);  z = bb(:,3);
    % 8 vertices in canonical order
    V = [x(1) y(1) z(1);  x(2) y(1) z(1);  x(2) y(2) z(1);  x(1) y(2) z(1);
         x(1) y(1) z(2);  x(2) y(1) z(2);  x(2) y(2) z(2);  x(1) y(2) z(2)];
    edges = [1 2; 2 3; 3 4; 4 1;     % bottom face
             5 6; 6 7; 7 8; 8 5;     % top face
             1 5; 2 6; 3 7; 4 8];    % verticals
    for e = 1:size(edges, 1)
        plot3(ax, V(edges(e,:), 1), V(edges(e,:), 2), V(edges(e,:), 3), ...
            '-', 'Color', [col alpha], 'LineWidth', lw, ...
            'HandleVisibility', 'off');
    end
end
