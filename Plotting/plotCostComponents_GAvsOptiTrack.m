function T = plotCostComponents_GAvsOptiTrack(varargin)
%PLOTCOSTCOMPONENTS_GAVSOPTITRACK  What each objective contributes to J.
%
%   T = plotCostComponents_GAvsOptiTrack('Name', Value, ...)
%
%   Stacked bars of the two weighted CF3 terms — resolution uncertainty
%   (Junc) and dynamic occlusion (Jocc) — for the Optimised GA Rig configuration and
%   for the Manually Posed Rig, across all four (target type x grid mode)
%   cells.
%
%   Two panels rather than one, because after the utopia/nadir
%   normalisation the GA sits near 0.05 and the ad-hoc rig above 1.5:
%   plotted on a shared linear axis the GA bars would be invisible, and a
%   log axis cannot carry a stack honestly. Each panel therefore has its
%   own scale, and the ad-hoc panel is annotated with the GA:ad-hoc ratio
%   so the comparison is still readable off the figure.
%
%   The dashed line at J = 1 is the nadir reference: both terms are scaled
%   so utopia -> 0 and nadir -> 1, and the weights sum to 1, so J = 1 is
%   the cost of a configuration that is at nadir in both objectives.
%
%   This figure exists because the corrected normalisation changed the
%   story. Under the old fixed constants occlusion accounted for almost
%   all of J and resolution uncertainty was numerically invisible; with
%   utopia/nadir scaling the two terms are comparable, so the split is
%   worth showing directly.
%
%   Name-Value parameters:
%     'NumCameras'  camera count. Default 7.
%     'Spacing'     grid spacing [m]. Default 1.0.
%     'LogFile'     master log. Default Results/Logs/GGA_RunsLog.mat.
%     'Breakdown'   a table from reportCostBreakdown, to skip recomputing.
%     'SaveAs'      output name without extension. Default auto.
%
%   Returns the breakdown table it plotted.
%
%   See also reportCostBreakdown, cf3Terms, plotCostField_GAvsOptiTrack.

    projectRoot = addProjectPaths();

    p = inputParser;
    addParameter(p, 'NumCameras', 7,   @isnumeric);
    addParameter(p, 'Spacing',    1.0, @isnumeric);
    addParameter(p, 'LogFile',    fullfile(projectRoot,'Results','Logs','GGA_RunsLog.mat'), @ischar);
    addParameter(p, 'Breakdown',  [],  @(x) isempty(x) || istable(x));
    addParameter(p, 'SaveAs',     '',  @ischar);
    parse(p, varargin{:});
    opts = p.Results;

    sty = gaPlotStyle();

    T = opts.Breakdown;
    if isempty(T)
        T = reportCostBreakdown('NumCameras', opts.NumCameras, ...
                                'Spacing', opts.Spacing, 'LogFile', opts.LogFile);
    end
    if isempty(T)
        error('plotCostComponents_GAvsOptiTrack:noData', ...
              'reportCostBreakdown returned no rows.');
    end

    % Two-line tick labels: the tick interpreter is tex, so \newline is
    % what actually breaks the line here (a literal newline is not).
    gmShort = replace(replace(string(T.GridMode), "Uniform", "Unif."), "Normal", "Norm.");
    labels  = strcat(string(T.TargetType), '\newline', gmShort);
    n      = height(T);

    uncCol = sty.CostFuncColors(1,:);   % resolution uncertainty
    occCol = sty.CostFuncColors(2,:);   % dynamic occlusion

    fig = figure('Name', 'CF3 cost composition', ...
        'Units', 'inches', ...
        'Position', [0.5, 0.5, sty.FigWidthFull * 1.25, sty.FigHeight + 0.9], ...
        'PaperPositionMode', 'auto', ...
        'Color', sty.BackgroundColor);

    tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    %% --- Panel 1: Optimised GA Rig ------------------------------------------------
    ax1 = nexttile(tl);
    b1 = drawStack(ax1, labels, [T.GA_Junc, T.GA_Jocc], uncCol, occCol, sty);
    ylabel(ax1, 'Weighted cost contribution', ...
        'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
    title(ax1, sprintf('Optimised GA Rig (k = %d)', opts.NumCameras), ...
        'FontWeight', 'normal', 'FontSize', sty.FontSizeTitle, 'FontName', sty.FontName);
    for i = 1:n
        text(ax1, i, T.GA_Total(i), sprintf(' %.3f', T.GA_Total(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', sty.FontSizeAnnot, 'FontName', sty.FontName);
    end
    ylim(ax1, [0, max(T.GA_Total) * 1.20]);

    %% --- Panel 2: ad-hoc rig ---------------------------------------------
    ax2 = nexttile(tl);
    drawStack(ax2, labels, [T.Adhoc_Junc, T.Adhoc_Jocc], uncCol, occCol, sty);
    title(ax2, 'Manually Posed Rig', ...
        'FontWeight', 'normal', 'FontSize', sty.FontSizeTitle, 'FontName', sty.FontName);
    hNadir = yline(ax2, 1, '--', 'nadir  (J = 1)', ...
        'Color', [0.30 0.30 0.30], 'LineWidth', 1.2, ...
        'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', ...
        'FontSize', sty.FontSizeAnnot, 'FontName', sty.FontName, ...
        'HandleVisibility', 'on', ...
        'DisplayName', 'nadir (J = 1)');
    for i = 1:n
        text(ax2, i, T.Adhoc_Total(i), ...
            sprintf(' %.2f\n(%.0f\\times)', T.Adhoc_Total(i), T.Ratio_Total(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', sty.FontSizeAnnot, 'FontName', sty.FontName, ...
            'Interpreter', 'tex');
    end
    ylim(ax2, [0, max(T.Adhoc_Total) * 1.28]);

    % One legend for the whole layout, parked below both panels, so it
    % cannot sit on top of a bar in either.
    lgd = legend([b1(1), b1(2), hNadir], ...
        {'J_{unc}  resolution uncertainty', ...
         'J_{occ}  dynamic occlusion', ...
         'nadir (J = 1)'}, ...
        'FontSize', sty.FontSizeLegend, 'Box', 'off', 'Orientation', 'horizontal');
    lgd.Layout.Tile = 'south';

    thesisTitle(tl, sprintf(['CF3 composition after utopia/nadir ' ...
        'normalisation (grid %.2f m, equal weights)  --  NOTE: independent ' ...
        'y-axis scales'], opts.Spacing), sty, 'MaxChars', 72);

    applyThesisStyle(fig);

    %% --- Export -----------------------------------------------------------
    if isempty(opts.SaveAs)
        outName = fullfile(projectRoot, 'figures', ...
            sprintf('CostComponents_GAvsOptiTrack_%dC_sp%.0fcm', ...
                    opts.NumCameras, opts.Spacing*100));
    else
        outName = opts.SaveAs;
    end
    exportThesisFigure(fig, outName, 'Background', sty.ExportBgColor, 'Quiet', true);
    fprintf('Saved: %s.{pdf,png}\n', outName);
end


%% ---- Local helper ------------------------------------------------------

function b = drawStack(ax, labels, vals, uncCol, occCol, sty)
    b = bar(ax, 1:numel(labels), vals, 0.6, 'stacked');
    b(1).FaceColor = uncCol;  b(1).EdgeColor = 'k';  b(1).LineWidth = 0.5;
    b(2).FaceColor = occCol;  b(2).EdgeColor = 'k';  b(2).LineWidth = 0.5;
    % Rotation off: MATLAB auto-rotates tick labels it thinks are too wide,
    % and a rotated two-line label overlaps its neighbour.
    set(ax, 'XTick', 1:numel(labels), 'XTickLabel', labels, ...
        'XTickLabelRotation', 0, ...
        'FontSize', sty.FontSizeTick, 'FontName', sty.FontName, ...
        'Box', 'on', 'TickDir', 'out');
    xlim(ax, [0.4, numel(labels) + 0.6]);
    grid(ax, 'on');
    ax.YGrid = 'on';
    ax.XGrid = 'off';
end
