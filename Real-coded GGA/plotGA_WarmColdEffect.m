function plotGA_WarmColdEffect(varargin)
% plotGA_WarmColdEffect  Compare warm-start vs cold-start performance.
%   plotGA_WarmColdEffect('Name', Value, ...)
%
%   For each camera count (and optionally split by target type), shows a
%   grouped bar chart of mean best cost for warm vs cold starts, with
%   individual data points overlaid and error bars (± 1 std).
%
%   Name-Value Parameters (passed through to loadGARuns, plus):
%     'SaveAs' - Output filename without extension (default: auto)
%     (all loadGARuns parameters are also accepted)

    p = inputParser;
    p.KeepUnmatched = true;
    addParameter(p, 'SaveAs', '', @ischar);
    parse(p, varargin{:});

    opts = p.Results;
    loaderArgs = [fieldnames(p.Unmatched), struct2cell(p.Unmatched)]';
    loaderArgs = loaderArgs(:)';

    [runs, filterDesc] = loadGARuns(loaderArgs{:});
    sty = gaPlotStyle();

    %% Check that both warm and cold runs exist
    warmFlags = [runs.WarmStart];
    if ~any(warmFlags) || ~all(~warmFlags | warmFlags)
        warning('Need both warm-start and cold-start runs for comparison.');
        if all(warmFlags)
            fprintf('  All %d runs are warm-start.\n', length(runs));
        else
            fprintf('  All %d runs are cold-start.\n', length(runs));
        end
        return;
    end

    uniqueCams = sort(unique([runs.NumCameras]));
    nCams = length(uniqueCams);

    %% Build summary table
    %  Rows: camera counts, Columns: [cold_mean, cold_std, warm_mean, warm_std]
    summaryData = zeros(nCams, 4);
    coldPts = cell(nCams, 1);
    warmPts = cell(nCams, 1);

    for c = 1:nCams
        cam = uniqueCams(c);
        coldMask = ([runs.NumCameras] == cam) & (~warmFlags);
        warmMask = ([runs.NumCameras] == cam) & (warmFlags);

        coldCosts = [runs(coldMask).BestCost];
        warmCosts = [runs(warmMask).BestCost];

        coldPts{c} = coldCosts;
        warmPts{c} = warmCosts;

        summaryData(c,:) = [mean(coldCosts), std(coldCosts), ...
                            mean(warmCosts), std(warmCosts)];
    end

    %% Plot
    fig = figure('Units', 'inches', ...
        'Position', [1 1 sty.FigWidthFull sty.FigHeight], ...
        'PaperPositionMode', 'auto', 'Color', 'w');
    ax = axes(fig);
    hold(ax, 'on');

    barWidth = 0.35;
    xPos = 1:nCams;

    % Cold bars
    hCold = bar(ax, xPos - barWidth/2, summaryData(:,1), barWidth, ...
        'FaceColor', sty.ColdColor, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    % Warm bars
    hWarm = bar(ax, xPos + barWidth/2, summaryData(:,3), barWidth, ...
        'FaceColor', sty.WarmColor, 'FaceAlpha', 0.6, 'EdgeColor', 'none');

    % Error bars
    errorbar(ax, xPos - barWidth/2, summaryData(:,1), summaryData(:,2), ...
        'k.', 'LineWidth', 0.8, 'CapSize', 4);
    errorbar(ax, xPos + barWidth/2, summaryData(:,3), summaryData(:,4), ...
        'k.', 'LineWidth', 0.8, 'CapSize', 4);

    % Individual data points (jittered)
    jitter = 0.06;
    for c = 1:nCams
        nC = length(coldPts{c});
        nW = length(warmPts{c});
        xCold = xPos(c) - barWidth/2 + jitter*(rand(1,nC)-0.5);
        xWarm = xPos(c) + barWidth/2 + jitter*(rand(1,nW)-0.5);

        plot(ax, xCold, coldPts{c}, 'o', ...
            'MarkerSize', 3, 'MarkerEdgeColor', sty.ColdColor*0.6, ...
            'MarkerFaceColor', sty.ColdColor*0.6, 'HandleVisibility', 'off');
        plot(ax, xWarm, warmPts{c}, 'o', ...
            'MarkerSize', 3, 'MarkerEdgeColor', sty.WarmColor*0.6, ...
            'MarkerFaceColor', sty.WarmColor*0.6, 'HandleVisibility', 'off');
    end

    hold(ax, 'off');

    %% Formatting
    set(ax, 'XTick', xPos, 'XTickLabel', arrayfun(@(x) sprintf('%d', x), uniqueCams, 'Uni', false));
    xlabel(ax, 'Number of Cameras', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
    ylabel(ax, 'Best Cost', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
    set(ax, 'FontSize', sty.FontSizeTick, 'FontName', sty.FontName, ...
        'Box', 'on', 'TickDir', 'out');
    grid(ax, 'on');
    set(ax, 'GridAlpha', 0.15);
    legend([hCold, hWarm], {'Cold-Start', 'Warm-Start'}, ...
        'Location', 'best', 'FontSize', sty.FontSizeLegend);

    %% Print improvement percentages
    fprintf('Warm-start improvement over cold-start:\n');
    for c = 1:nCams
        coldMean = summaryData(c,1);
        warmMean = summaryData(c,3);
        pctImprove = 100 * (coldMean - warmMean) / coldMean;
        fprintf('  %dC: Cold=%.4f, Warm=%.4f (%.1f%% improvement)\n', ...
            uniqueCams(c), coldMean, warmMean, pctImprove);
    end

    %% Export
    if isempty(opts.SaveAs)
        outName = sprintf('GA_WarmCold_%s', filterDesc);
    else
        outName = opts.SaveAs;
    end
    exportgraphics(fig, [outName '.pdf'], 'ContentType', 'vector');
    fprintf('Saved: %s.pdf\n', outName);
end