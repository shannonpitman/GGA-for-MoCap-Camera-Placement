function paths = plotOrientationSensitivity(sweep, projectRoot, ts)
%PLOTORIENTATIONSENSITIVITY  Figures for the orientation sensitivity study.
%
%   plotOrientationSensitivity(sweep, projectRoot, ts) produces four
%   figures from the struct returned by orientationSensitivity:
%
%     1  cost vs orientation grid step, one panel per cost function
%     2  CF3 deviation from the unsnapped optimum, with tolerance band
%     3  how far each camera had to turn to reach the grid, split into the
%        roll component (free — it does not move the optical axis) and the
%        optical-axis component, which is the part that moves coverage.
%        Roll is reported modulo 180 degrees, since a 180-degree roll maps
%        the sensor rectangle onto itself.
%     4  Part A: per-camera roll before and after uprighting, and how close
%        the unconstrained optimum already sat to the 15-degree grid
%
%   Colour encodes the snapped gene set: all three angles, pointing only
%   (alpha, beta) or roll only (gamma). Only the GA configuration appears —
%   the ad-hoc rig is already installed and is never snapped, so plotting a
%   quantisation curve for it would be meaningless. Where a config axis
%   exists at all it is drawn with line style, so a sweep carrying more
%   than one configuration still plots sensibly.
%
%   Files land in figures/Sensitivity/Orientation/ as PDF + PNG via
%   exportThesisFigure, so they drop straight into the thesis figure set.
%
%   See also: orientationSensitivity, plotSpacingSensitivity.

    if nargin < 3 || isempty(ts)
        ts = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    end
    if nargin < 2 || isempty(projectRoot)
        projectRoot = addProjectPaths();
    end

    sty    = gaPlotStyle();
    figDir = fullfile(projectRoot, 'figures', 'Sensitivity', 'Orientation');
    if ~isfolder(figDir), mkdir(figDir); end

    cfgNames = unique({sweep.snap.config}, 'stable');
    nC       = numel(cfgNames);
    variants = sweep.variants;
    nV       = numel(variants);
    steps    = sweep.steps;

    % Colour = gene set. Line style = configuration (usually just one).
    vColors    = sty.CostFuncColors;
    if size(vColors,1) < nV, vColors = lines(nV); end
    vMarkers   = {'o','s','d'};
    lineStyles = {'-', '--', ':'};

    tag = sprintf('%s_GM%d_sp%03.0fcm_%dC', sweep.modeTag, sweep.gridMode, ...
                  sweep.spacing*100, sweep.numCams);
    paths = {};

    %% ---- Figure 1: cost vs grid step ---------------------------------
    cfFields = {'CF1', 'CF2', 'CF3'};
    cfLabels = {'Resolution uncertainty (raw)', ...
                'Dynamic occlusion (raw)', ...
                'Combined CF3 (normalised)'};

    fig1 = figure('Units','inches', 'Position',[1 1 sty.FigWidthDouble sty.FigHeight], 'Color','w');
    tl = tiledlayout(1, 3, 'TileSpacing','compact', 'Padding','compact');
    for f = 1:3
        nexttile; hold on;
        for c = 1:nC
            for v = 1:nV
                rows = snapRows(sweep, cfgNames{c}, variants(v).key);
                if isempty(rows), continue; end
                y = arrayfun(@(r) r.cost.(cfFields{f}), rows);
                plot([rows.stepDeg], y, [lineStyles{min(c,end)} vMarkers{min(v,end)}], ...
                     'LineWidth', sty.LineWidth, 'MarkerSize', sty.MarkerSize, ...
                     'Color', vColors(v,:), 'MarkerFaceColor', vColors(v,:), ...
                     'DisplayName', seriesName(cfgNames, c, variants(v).key));
            end
            ref = uprightValue(sweep, cfgNames{c}, 'upright-flip', cfFields{f});
            yline(ref, '-', 'Color', [0.45 0.45 0.45], 'LineWidth', sty.LineWidth, ...
                  'HandleVisibility', 'off');
        end
        xline(15, ':', 'Color', [0.3 0.3 0.3], 'LineWidth', sty.LineWidthThin, ...
              'HandleVisibility', 'off');
        set(gca, 'XScale', 'log', 'YScale', 'log', ...
                 'XTick', steps, 'XTickLabel', compose('%g', steps), ...
                 'XTickLabelRotation', 45);
        xlabel('Orientation grid step [deg]');
        ylabel('Cost');
        title(cfLabels{f}, 'FontWeight', 'normal', 'FontSize', sty.FontSizeTitle);
        grid on; box on;
    end
    lg = legend('Box','off', 'FontSize', sty.FontSizeLegend, 'NumColumns', 3);
    lg.Layout.Tile = 'south';
    title(tl, sprintf('Cost vs orientation quantisation — %s, %d cameras', ...
        sweep.modeTag, sweep.numCams), 'FontWeight','bold');
    subtitle(tl, 'Grey line = unsnapped upright optimum', ...
        'FontSize', sty.FontSizeAnnot);
    applyThesisStyle(fig1);
    paths{end+1} = exportThesisFigure(fig1, fullfile(figDir, ...
        sprintf('orientationsnap_cost_%s_%s', tag, ts)));

    %% ---- Figure 2: CF3 deviation with tolerance band -----------------
    fig2 = figure('Units','inches', 'Position',[1 1 sty.FigWidthFull sty.FigHeight], 'Color','w');
    hold on;
    tol = sweep.tolerancePct;
    xl  = [min(steps)*0.8, max(steps)*1.25];
    if ~isempty(tol)
        fill([xl fliplr(xl)], [-tol -tol tol tol], [0.85 0.9 0.85], ...
             'EdgeColor','none', 'FaceAlpha', 0.55, ...
             'DisplayName', sprintf('\\pm%.0f%% tolerance', tol));
    end
    for c = 1:nC
        ref = uprightValue(sweep, cfgNames{c}, 'upright-flip', 'CF3');
        for v = 1:nV
            rows = snapRows(sweep, cfgNames{c}, variants(v).key);
            if isempty(rows), continue; end
            y      = arrayfun(@(r) r.cost.CF3, rows);
            devPct = 100 * (y - ref) / abs(ref);
            plot([rows.stepDeg], devPct, [lineStyles{min(c,end)} vMarkers{min(v,end)}], ...
                 'LineWidth', sty.LineWidth, 'MarkerSize', sty.MarkerSize, ...
                 'Color', vColors(v,:), 'MarkerFaceColor', vColors(v,:), ...
                 'DisplayName', seriesName(cfgNames, c, variants(v).key));
        end
    end
    yline(0, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', sty.LineWidthThin, 'HandleVisibility','off');
    xline(15, ':', '15\circ', 'Color', [0.2 0.2 0.2], 'LineWidth', sty.LineWidth, ...
          'LabelVerticalAlignment','bottom', 'HandleVisibility','off');
    set(gca, 'XScale','log', 'XTick', steps, 'XTickLabel', compose('%g', steps), ...
             'XTickLabelRotation', 45);
    xlim(xl);
    xlabel('Orientation grid step [deg]');
    ylabel('CF3 change vs unsnapped optimum [%]');
    title(sprintf('Cost of mounting the optimised rig to a discrete angular grid — %s, %d cameras', ...
        sweep.modeTag, sweep.numCams), 'FontWeight','bold');
    legend('Location','northwest', 'Box','off', 'FontSize', sty.FontSizeLegend);
    grid on; box on;
    applyThesisStyle(fig2);
    paths{end+1} = exportThesisFigure(fig2, fullfile(figDir, ...
        sprintf('orientationsnap_deviation_%s_%s', tag, ts)));

    %% ---- Figure 3: how far the cameras had to turn -------------------
    fig3 = figure('Units','inches', 'Position',[1 1 sty.FigWidthDouble sty.FigHeight], 'Color','w');
    tl = tiledlayout(1, 2, 'TileSpacing','compact', 'Padding','compact');

    nexttile; hold on;
    for c = 1:nC
        for v = 1:nV
            rows = snapRows(sweep, cfgNames{c}, variants(v).key);
            if isempty(rows), continue; end
            maxR = arrayfun(@(r) max(abs(r.snapReport.RollShiftDeg), [], 'omitnan'), rows);
            plot([rows.stepDeg], maxR, [lineStyles{min(c,end)} vMarkers{min(v,end)}], ...
                 'LineWidth', sty.LineWidth, 'MarkerSize', sty.MarkerSize, ...
                 'Color', vColors(v,:), 'MarkerFaceColor', vColors(v,:), ...
                 'DisplayName', seriesName(cfgNames, c, variants(v).key));
        end
    end
    set(gca, 'XScale','log', 'XTick', steps, 'XTickLabel', compose('%g', steps), ...
             'XTickLabelRotation', 45);
    xlabel('Orientation grid step [deg]');
    ylabel('Max roll shift [deg]');
    title('Roll — free, does not move the axis', 'FontWeight','normal', ...
          'FontSize', sty.FontSizeTitle);
    grid on; box on;

    nexttile; hold on;
    for c = 1:nC
        for v = 1:nV
            rows = snapRows(sweep, cfgNames{c}, variants(v).key);
            if isempty(rows), continue; end
            maxA = arrayfun(@(r) max(r.snapReport.AxisShiftDeg), rows);
            plot([rows.stepDeg], maxA, [lineStyles{min(c,end)} vMarkers{min(v,end)}], ...
                 'LineWidth', sty.LineWidth, 'MarkerSize', sty.MarkerSize, ...
                 'Color', vColors(v,:), 'MarkerFaceColor', vColors(v,:), ...
                 'DisplayName', seriesName(cfgNames, c, variants(v).key));
        end
    end
    set(gca, 'XScale','log', 'XTick', steps, 'XTickLabel', compose('%g', steps), ...
             'XTickLabelRotation', 45);
    xlabel('Orientation grid step [deg]');
    ylabel('Max optical-axis shift [deg]');
    title('Optical axis — this drives coverage', 'FontWeight','normal', ...
          'FontSize', sty.FontSizeTitle);
    % 'all' is hidden exactly beneath 'pointing' here, and that coincidence
    % is the whole point: the optical axis is a function of alpha and beta
    % alone, so snapping gamma as well changes nothing on this axis.
    text(0.04, 0.94, ['''all'' lies exactly under ''pointing'':' newline ...
         'roll cannot move the optical axis'], ...
         'Units','normalized', 'HorizontalAlignment','left', ...
         'VerticalAlignment','top', 'Color', [0.3 0.3 0.3], ...
         'FontSize', sty.FontSizeAnnot);
    grid on; box on;

    lg = legend('Box','off', 'FontSize', sty.FontSizeLegend, 'NumColumns', 3);
    lg.Layout.Tile = 'south';
    title(tl, sprintf('Displacement imposed by snapping — %s, %d cameras', ...
        sweep.modeTag, sweep.numCams), 'FontWeight','bold');
    applyThesisStyle(fig3);
    paths{end+1} = exportThesisFigure(fig3, fullfile(figDir, ...
        sprintf('orientationsnap_displacement_%s_%s', tag, ts)));

    %% ---- Figure 4: roll canonicalisation + grid proximity ------------
    fig4 = figure('Units','inches', 'Position',[1 1 sty.FigWidthFull sty.FigHeightTall], 'Color','w');
    tl = tiledlayout(2, 1, 'TileSpacing','compact', 'Padding','compact');

    barColor = sty.CameraColors(1,:);
    nexttile; hold on;
    camIdx = 1:sweep.numCams;
    before = uprightInfo(sweep, cfgNames{1}, 'baseline');
    after  = uprightInfo(sweep, cfgNames{1}, 'upright-flip');
    bar(camIdx, before.TiltDeg, 0.55, ...
        'FaceColor', barColor, 'FaceAlpha', 0.35, 'EdgeColor', barColor, ...
        'DisplayName', 'As optimised');
    plot(camIdx, after.TiltDeg, 'o', ...
         'MarkerSize', sty.MarkerSizeLg, 'MarkerFaceColor', barColor, ...
         'MarkerEdgeColor', 'k', 'LineStyle','none', ...
         'DisplayName', 'After 180\circ flip');
    yline( 90, '--', 'Color', [0.75 0.2 0.2], 'LineWidth', sty.LineWidthThin, 'HandleVisibility','off');
    yline(-90, '--', 'Color', [0.75 0.2 0.2], 'LineWidth', sty.LineWidthThin, 'HandleVisibility','off');
    % Label on the left, where no camera bar reaches above +90.
    text(0.5, 130, ' inverted', 'Color', [0.75 0.2 0.2], ...
         'FontSize', sty.FontSizeAnnot, 'HorizontalAlignment','left');
    ylim([-195 195]); yticks(-180:90:180);
    xlim([0.4, sweep.numCams + 0.6]); xticks(camIdx);
    xlabel('Camera'); ylabel('Roll about optical axis [deg]');
    title('Part A — roll relative to level, before and after uprighting', ...
          'FontWeight','normal', 'FontSize', sty.FontSizeTitle);
    % Bars run well above the axis but never below -90, so the bottom of
    % the panel is the one region guaranteed clear of data.
    legend('Location','southeast', 'Box','off', 'FontSize', sty.FontSizeLegend, ...
           'NumColumns', 2);
    grid on; box on;

    nexttile; hold on;
    rows = snapRows(sweep, cfgNames{1}, 'all');
    if isempty(rows), rows = snapRows(sweep, cfgNames{1}, variants(1).key); end
    at15 = rows([rows.stepDeg] == 15);
    if ~isempty(at15)
        r = at15(1).snapReport.GridResidualDeg(:);
        histogram(r, 0:0.75:7.5, 'FaceColor', barColor, 'FaceAlpha', 0.55, ...
                  'EdgeColor', barColor, 'HandleVisibility','off');
        xline(7.5/2, '--', 'uniform mean', 'Color', [0.4 0.4 0.4], ...
              'LineWidth', sty.LineWidthThin, 'HandleVisibility','off');
    end
    xlabel('Distance from nearest 15\circ multiple [deg]');
    ylabel('Euler genes');
    title('Part B — how close the unconstrained optimum already sat to the grid', ...
          'FontWeight','normal', 'FontSize', sty.FontSizeTitle);
    grid on; box on;

    title(tl, sprintf('Camera roll and angular quantisation — %s, %d cameras', ...
        sweep.modeTag, sweep.numCams), 'FontWeight','bold');
    applyThesisStyle(fig4);
    paths{end+1} = exportThesisFigure(fig4, fullfile(figDir, ...
        sprintf('orientation_uprighting_%s_%s', tag, ts)));
end


%% ======================================================================
%  LOCAL HELPERS
%  ======================================================================

function name = seriesName(cfgNames, c, variantKey)
% With the usual single configuration the config name adds nothing to every
% legend entry, so it is only prefixed when there is more than one.
    if numel(cfgNames) == 1
        name = variantKey;
    else
        name = sprintf('%s — %s', cfgNames{c}, variantKey);
    end
end

function rows = snapRows(sweep, configName, variantKey)
    mask = strcmp({sweep.snap.config}, configName);
    if nargin >= 3 && ~isempty(variantKey)
        mask = mask & strcmp({sweep.snap.variant}, variantKey);
    end
    rows = sweep.snap(mask);
    if isempty(rows), return; end
    [~, order] = sort([rows.stepDeg]);
    rows = rows(order);
end

function v = uprightValue(sweep, configName, variant, field)
    row = sweep.upright(strcmp({sweep.upright.config}, configName) & ...
                        strcmp({sweep.upright.variant}, variant));
    v = row(1).cost.(field);
end

function info = uprightInfo(sweep, configName, variant)
    row = sweep.upright(strcmp({sweep.upright.config}, configName) & ...
                        strcmp({sweep.upright.variant}, variant));
    info = row(1).info;
end
