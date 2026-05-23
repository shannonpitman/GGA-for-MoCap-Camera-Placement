function plotSpacingSensitivity(sweep, projectRoot, ts, tolerancePct, recs)
%PLOTSPACINGSENSITIVITY  Three diagnostic figures from a sweep struct.
%
%   plotSpacingSensitivity(sweep, projectRoot, ts) generates:
%     1. cost vs spacing (one panel per CF)
%     2. relative deviation from finest grid (% change)
%     3. evaluation time + grid-point count vs spacing
%
%   plotSpacingSensitivity(..., tolerancePct, recs) additionally overlays:
%     - tolerance band (+/- tolerancePct) on the deviation plot
%     - vertical line at the recommended spacing for each (config, CF)
%
%   Outputs are saved as PNGs into figures/Sensitivity/<modeTag>/ where
%   modeTag is taken from sweep.modeTag (e.g. 'UAV', 'UGV').

    if nargin < 4 || isempty(tolerancePct), tolerancePct = []; end
    if nargin < 5, recs = []; end

    modeTag = '';
    if isfield(sweep, 'modeTag') && ~isempty(sweep.modeTag)
        modeTag = sweep.modeTag;
    end

    figDir = fullfile(projectRoot, 'figures', 'Sensitivity', modeTag);
    if ~isfolder(figDir), mkdir(figDir); end

    results     = sweep.results;
    cfNames     = sweep.cfNames;
    cfLabels    = sweep.cfLabels;
    configs     = sweep.configs;
    configNames = {configs.name};
    nC = numel(configNames);
    nF = numel(cfNames);

    colors = lines(nC);
    markers = {'o','s','d','^'};

    %% Fig 1 — Cost vs spacing
    fig1 = figure('Position', [80 80 1500 440], 'Color', 'w');
    tl = tiledlayout(1, nF, 'TileSpacing', 'compact', 'Padding', 'compact');
    for f = 1:nF
        nexttile; hold on;
        for c = 1:nC
            sub = filterRows(results, configNames{c}, cfNames{f});
            plot(sub.spacing, sub.cost, ['-' markers{min(c,end)}], ...
                 'LineWidth', 1.6, 'MarkerSize', 6, ...
                 'Color', colors(c,:), 'MarkerFaceColor', colors(c,:), ...
                 'DisplayName', configNames{c});
        end
        xlabel('Grid spacing [m]');
        ylabel('Cost');
        title(cleanTitle(cfLabels{f}), 'FontWeight', 'normal');
        set(gca, 'XScale', 'log');
        grid on; box on;
        if f == nF
            legend('Location', 'best', 'Box', 'off');
        end
    end
    titleStr = sprintf('Cost vs grid spacing: %s, volume %g x %g x %g m', ...
        upper(modeTag), diff(sweep.volume(1,:)), diff(sweep.volume(2,:)), diff(sweep.volume(3,:)));
    title(tl, titleStr, 'FontWeight', 'bold');
    applyThesisStyle(fig1);
    f1Path = fullfile(figDir, sprintf('cost_vs_spacing_%s.png', ts));
    exportgraphics(fig1, f1Path, 'Resolution', 150, 'BackgroundColor', 'white');

    %% Fig 2 — Relative deviation from finest-grid reference
    fig2 = figure('Position', [80 80 1500 440], 'Color', 'w');
    tl = tiledlayout(1, nF, 'TileSpacing', 'compact', 'Padding', 'compact');
    for f = 1:nF
        nexttile; hold on;
        for c = 1:nC
            sub = filterRows(results, configNames{c}, cfNames{f});
            ref = sub.cost(1);
            relPct = (sub.cost - ref) ./ abs(ref) * 100;
            plot(sub.spacing, relPct, ['-' markers{min(c,end)}], ...
                 'LineWidth', 1.6, 'MarkerSize', 6, ...
                 'Color', colors(c,:), 'MarkerFaceColor', colors(c,:), ...
                 'DisplayName', configNames{c});
        end

        % Tolerance band
        if ~isempty(tolerancePct)
            yline(0, ':k');
            yline( tolerancePct, '--', sprintf('+%.1f%%', tolerancePct), ...
                'LabelHorizontalAlignment','left', 'Color',[0.4 0.4 0.4]);
            yline(-tolerancePct, '--', sprintf('-%.1f%%', tolerancePct), ...
                'LabelHorizontalAlignment','left', 'Color',[0.4 0.4 0.4]);
        end

        % Recommended-spacing markers
        if ~isempty(recs)
            for c = 1:nC
                idx = strcmp({recs.config}, configNames{c}) & ...
                      strcmp({recs.cf},     cfNames{f});
                if any(idx)
                    xrec = recs(idx).recommendedSpacing;
                    xline(xrec, ':', sprintf('rec %.2gm', xrec), ...
                          'Color', colors(c,:), 'LabelOrientation','horizontal', ...
                          'LabelVerticalAlignment','top', 'LineWidth', 1.0);
                end
            end
        end

        xlabel('Grid spacing [m]');
        ylabel('Deviation from finest grid [%]');
        title(cleanTitle(cfLabels{f}), 'FontWeight', 'normal');
        set(gca, 'XScale', 'log');
        grid on; box on;
        if f == nF
            legend('Location', 'best', 'Box', 'off');
        end
    end
    devTitle = sprintf('Relative deviation from finest grid (%g m): %s', ...
        sweep.spacings(1), upper(modeTag));
    title(tl, devTitle, 'FontWeight', 'bold');
    applyThesisStyle(fig2);
    f2Path = fullfile(figDir, sprintf('relative_deviation_%s.png', ts));
    exportgraphics(fig2, f2Path, 'Resolution', 150, 'BackgroundColor', 'white');

    %% Fig 3 — Evaluation cost
    fig3 = figure('Position', [80 80 900 440], 'Color', 'w');
    yyaxis left
    hold on;
    for c = 1:nC
        sub = filterRows(results, configNames{c}, 'CF3_combined');
        plot(sub.spacing, sub.evalTime, ['-' markers{min(c,end)}], ...
             'LineWidth', 1.6, 'MarkerSize', 6, ...
             'Color', colors(c,:), 'MarkerFaceColor', colors(c,:), ...
             'DisplayName', sprintf('%s: eval time', configNames{c}));
    end
    set(gca, 'XScale', 'log', 'YScale', 'log');
    ylabel('CF3 evaluation time [s]');

    yyaxis right
    sub = filterRows(results, configNames{1}, 'CF3_combined');
    plot(sub.spacing, sub.numPoints, ':k', 'LineWidth', 1.2, ...
         'DisplayName', 'Grid points');
    set(gca, 'YScale', 'log');
    ylabel('# target points');

    xlabel('Grid spacing [m]');
    title(sprintf('Evaluation cost vs grid spacing: %s (CF3 = uncert + occl)', upper(modeTag)));
    grid on; box on;
    legend('Location', 'best', 'Box', 'off');

    applyThesisStyle(fig3);
    f3Path = fullfile(figDir, sprintf('eval_time_%s.png', ts));
    exportgraphics(fig3, f3Path, 'Resolution', 150, 'BackgroundColor', 'white');

    fprintf('Figures saved to: %s\n', figDir);
end


function s = cleanTitle(s)
% Replace em-dash / en-dash with a colon for thesis-figure consistency.
    s = strrep(s, ' — ', ': ');
    s = strrep(s, ' – ', ': ');
    s = strrep(s, '—', ':');
    s = strrep(s, '–', ':');
end


function sub = filterRows(results, configName, cfn)
    rows = results( strcmp({results.config}, configName) & ...
                    strcmp({results.cf},     cfn) );
    [sp, ord] = sort([rows.spacing]);
    sub.spacing   = sp(:);
    sub.cost      = reshape([rows(ord).cost], [], 1);
    sub.evalTime  = reshape([rows(ord).evalTime], [], 1);
    sub.numPoints = reshape([rows(ord).numPoints], [], 1);
end
