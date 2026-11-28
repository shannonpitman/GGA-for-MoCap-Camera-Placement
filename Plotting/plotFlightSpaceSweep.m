function R = plotFlightSpaceSweep(varargin)
%PLOTFLIGHTSPACESWEEP  OpenMV vs OptiTrack on tripods in the measured flight space.
%
%   R = plotFlightSpaceSweep()                 % layouts: full-coverage rule
%   R = plotFlightSpaceSweep('LayoutCF', 3)    % layouts: best CF3 run
%   R = plotFlightSpaceSweep('Presets', {'openmv_cab_2m', 'openmv_floor_2m', 'openmv_2m'}, ...
%                            'Tag', 'cabinets')
%
%   Uses the 'openmv_2m' (6-12 cameras) and 'optitrack_tripod' (6-7)
%   sweeps: 5 x 6 m flight space, targets 0.2-2 m, tripods 0.75-1.75 m in
%   the 1.5 m band around it. Every run's best rig is re-evaluated for its
%   raw resolution uncertainty [mm] and dynamic occlusion angle [deg], so
%   the comparison is free of the per-instance CF3 normalisation.
%
%   Figures (figures/flightspace/, PDF + PNG):
%     flightspace_sweep            raw uncertainty and occlusion against
%                                  camera count; CF2 (occlusion-only) and
%                                  CF3 (combined) runs, median and interquartile range;
%                                  OptiTrack on tripods as reference
%     flightspace_layouts_<sel>    top view of the chosen rig per camera count
%                                  (arrows = optical axis, colour = height)
%     flightspace_occlusion_<sel>  per-point occlusion angle of those rigs,
%                                  averaged over the three height layers
%
%   Name-value:
%     'LayoutCF'  which rig is drawn per camera count:
%                 0 (default, <sel> = 'fullcov') the lowest-occlusion rig
%                   among CF2 and CF3 runs that see every target point
%                   with at least two cameras. CF2 runs can leave a few
%                   single-camera points (0.75 m uncertainty penalty each),
%                   and CF3 at 6-7 cameras trades large occlusion for tiny
%                   resolution gains, so neither best-cost rig alone is the
%                   best rig in raw terms.
%                 1, 2, 3 (<sel> = 'CF<k>') the lowest-cost run of that CF
%     'Presets'   runConfig presets to compare (default openmv_2m and
%                 optitrack_tripod); the first is drawn with filled markers
%     'Names'     display names, one per preset (default: preset names)
%     'Tag'       figure file prefix (default 'flightspace')
%     'OutDir'    default <root>/figures/flightspace
%
%   R is a table: Preset, NumCameras, CostFunction, RunFilename, BestCost,
%   RawUnc_mm, RawOcc_deg, TwoPlusPct.

    root = addProjectPaths();
    p = inputParser;
    addParameter(p, 'LayoutCF', 0, @(x) any(x == [0 1 2 3]));
    addParameter(p, 'OutDir', fullfile(root, 'figures', 'flightspace'), @ischar);
    addParameter(p, 'Presets', {'openmv_2m', 'optitrack_tripod'}, @iscellstr);
    addParameter(p, 'Names', {}, @iscellstr);
    addParameter(p, 'Tag', 'flightspace', @ischar);
    parse(p, varargin{:});
    o = p.Results;
    if ~isfolder(o.OutDir), mkdir(o.OutDir); end
    s = gaPlotStyle();

    presets = o.Presets;  nPre = numel(presets);
    if ~isempty(o.Names)
        names = o.Names;
    elseif isequal(presets, {'openmv_2m', 'optitrack_tripod'})
        names = {'OpenMV RT1062', 'OptiTrack'};
    else
        names = strrep(presets, '_', '\_');
    end
    cfgs = cellfun(@runConfig, presets, 'UniformOutput', false);

    %% Raw components of every run's best rig
    rows = {};
    for k = 1:numel(presets)
        L = loadGARuns('Preset', presets{k});
        for i = 1:numel(L)
            sd = load(resolveRunPath(L(i).RunFilename, L(i).NumCameras), 'saveData').saveData;
            ch = runChromosome(sd);  sp = sd.Specifications;
            b = costBreakdown(ch, sp);
            [~, st] = perTargetCoverage(ch, sp);
            rows(end+1, :) = {presets{k}, L(i).NumCameras, L(i).CostFunctionType, ...
                L(i).RunFilename, L(i).BestCost, 1000*b.RawUnc, b.RawOcc, st.twoPlusPct}; %#ok<AGROW>
        end
    end
    R = cell2table(rows, 'VariableNames', {'Preset', 'NumCameras', 'CostFunction', ...
        'RunFilename', 'BestCost', 'RawUnc_mm', 'RawOcc_deg', 'TwoPlusPct'});

    %% Figure 1: camera-count sweep
    fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 s.FigWidthFull s.FigHeightTall]);
    tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    metrics = {'RawOcc_deg', 'RawUnc_mm'};
    ylabs   = {'Dynamic occlusion [deg]', 'Resolution uncertainty [mm]'};
    cfs = [2 3];  mk = {'o', 's', '^', 'd', 'v'};
    off = linspace(-0.25, 0.25, 2*nPre);
    allN = unique(R.NumCameras);
    for m = 1:2
        ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on'); h = gobjects(0); lab = {};
        for c = 1:2
            col = s.CostFuncColors(cfs(c), :);
            for k = 1:nPre
                T = R(strcmp(R.Preset, presets{k}) & R.CostFunction == cfs(c), :);
                if isempty(T), continue; end
                n = unique(T.NumCameras).';
                med = arrayfun(@(x) median(T.(metrics{m})(T.NumCameras == x)), n);
                lo  = arrayfun(@(x) prctile(T.(metrics{m})(T.NumCameras == x), 25), n);
                hi  = arrayfun(@(x) prctile(T.(metrics{m})(T.NumCameras == x), 75), n);
                x = n + off((k - 1)*2 + c);
                if k == 1
                    h(end+1) = errorbar(ax, x, med, med - lo, hi - med, ['-' mk{k}], 'Color', col, ...
                        'MarkerFaceColor', col, 'LineWidth', s.LineWidth, 'MarkerSize', s.MarkerSize); %#ok<AGROW>
                else
                    ls = mk{k};  if nPre > 2, ls = [':' mk{k}]; end
                    h(end+1) = errorbar(ax, x, med, med - lo, hi - med, ls, 'Color', col, ...
                        'MarkerFaceColor', 'w', 'LineWidth', s.LineWidth, 'MarkerSize', s.MarkerSizeLg); %#ok<AGROW>
                end
                lab{end+1} = sprintf('%s, %s', names{k}, s.CostFuncShort{cfs(c)}); %#ok<AGROW>
            end
        end
        ylabel(ax, ylabs{m}); xlim(ax, [min(allN) - 0.5, max(allN) + 0.5]); xticks(ax, allN.');
        ylim(ax, [0, ax.YLim(2)]);
        if m == 1, legend(ax, h, lab, 'Location', 'northeast', 'FontSize', s.FontSizeLegend); end
        if m == 2, xlabel(ax, 'Number of cameras'); end
    end
    title(tl, 'Tripod rigs, 5 \times 6 m flight space (median, IQR)', ...
        'FontSize', s.FontSizeTitle);
    applyThesisStyle(fig);
    exportThesisFigure(fig, fullfile(o.OutDir, [o.Tag '_sweep']));

    %% Best rig per camera count for the layout figures
    sets = {};
    for k = 1:nPre
        cfg = cfgs{k};
        for nc = cfg.CameraRange
            if o.LayoutCF == 0
                T = R(strcmp(R.Preset, presets{k}) & R.NumCameras == nc & ...
                      ismember(R.CostFunction, [2 3]) & R.TwoPlusPct >= 100 - 1e-9, :);
                if isempty(T), continue; end
                [~, i] = min(T.RawOcc_deg);
            else
                T = R(strcmp(R.Preset, presets{k}) & R.NumCameras == nc & R.CostFunction == o.LayoutCF, :);
                if isempty(T), continue; end
                [~, i] = min(T.BestCost);
            end
            sd = load(resolveRunPath(T.RunFilename{i}, nc), 'saveData').saveData;
            sets(end+1, :) = {sprintf('%d %s', nc, names{k}), runChromosome(sd), sd.Specifications, ...
                T.RawUnc_mm(i), T.RawOcc_deg(i), T.TwoPlusPct(i), T.CostFunction(i), ...
                k, find(allN == nc)}; %#ok<AGROW>
        end
    end
    if o.LayoutCF == 0
        sel = 'fullcov';  selTitle = 'Lowest-occlusion full-coverage rigs';
    else
        sel = sprintf('CF%d', o.LayoutCF);  selTitle = sprintf('Best %s rigs', s.CostFuncShort{o.LayoutCF});
    end
    % Grid: presets side by side; with one or two presets several camera
    % counts share a row (up to 4 columns) so the tiles stay readable
    nP = size(sets, 1);
    perRow = max(1, floor(4 / nPre));            % camera counts per row
    nCol = perRow * nPre;  nRow = ceil(numel(allN) / perRow);
    tile = @(j) (ceil(sets{j, 9} / perRow) - 1)*nCol + mod(sets{j, 9} - 1, perRow)*nPre + sets{j, 8};
    rowH = (s.FigWidthFull - 0.9) / nCol * 1.15 + 0.75;   % equal-aspect tile + title lines
    xyLim = max(cell2mat(cellfun(@(c) c.CamUpperBounds(1:2), cfgs(:), 'UniformOutput', false)), [], 1);
    tz = [inf -inf];
    for k = 1:nPre
        M = cfgs{k}.Mount;  hz = M.TripodHeight;
        if isfield(M, 'EndHeight'), hz = [hz M.EndHeight]; end %#ok<AGROW>
        if isfield(M, 'Cabinet'), hz = [hz M.Cabinet.Height + M.TripodHeight]; end %#ok<AGROW>
        tz = [min(tz(1), min(hz)) max(tz(2), max(hz))];
    end

    %% Figure 2: layouts
    fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 s.FigWidthFull rowH*nRow]);
    tl = tiledlayout(fig, nRow, nCol, 'TileSpacing', 'compact', 'Padding', 'compact');
    for j = 1:nP
        ax = nexttile(tl, tile(j)); hold(ax, 'on');
        cfg = cfgs{sets{j, 8}};
        fx = cfg.Volume(1, :);  fy = cfg.Volume(2, :);
        bx = [-1 1] * min(abs(cfg.CamLowerBounds(1)), 4);  by = [-1 1] * abs(cfg.CamLowerBounds(2));
        rectangle(ax, 'Position', [bx(1) by(1) diff(bx) diff(by)], 'FaceColor', [0.93 0.93 0.93], 'EdgeColor', [0.6 0.6 0.6]);
        if isfield(cfg.Mount, 'Cabinet')
            cb = cfg.Mount.Cabinet;
            for sx = [-1 1]
                xs = sort(sx * cb.XDepth);
                rectangle(ax, 'Position', [xs(1) cb.YSpan(1) diff(xs) diff(cb.YSpan)], ...
                    'FaceColor', [0.82 0.74 0.62], 'EdgeColor', [0.5 0.4 0.3]);
            end
        end
        rectangle(ax, 'Position', [fx(1) fy(1) diff(fx) diff(fy)], 'FaceColor', 'w', 'EdgeColor', 'k');
        P = reshape(sets{j, 2}, 6, []).';
        for c = 1:size(P, 1)
            Rm = genesToRotm(P(c, 4:6));  d = Rm(:, 3);  dxy = d(1:2) / max(norm(d(1:2)), eps);
            quiver(ax, P(c, 1), P(c, 2), 1.1*dxy(1), 1.1*dxy(2), 0, 'Color', [0.3 0.3 0.3], ...
                'LineWidth', 1, 'MaxHeadSize', 0.6);
        end
        scatter(ax, P(:, 1), P(:, 2), 36, P(:, 3), 'filled', 'MarkerEdgeColor', 'k');
        clim(ax, tz); colormap(ax, parula);
        axis(ax, 'equal'); xlim(ax, xyLim(1) * [-1 1] + [-0.3 0.3]); ylim(ax, xyLim(2) * [-1 1] + [-0.3 0.3]);
        ax.XTick = -4:2:4; ax.YTick = -4:2:4; box(ax, 'on');
        title(ax, sprintf('%s', sets{j, 1}), 'FontSize', s.FontSizeAnnot);
        st = {sprintf('%.1f mm, %.0f deg', sets{j, 4}, sets{j, 5})};
        if sets{j, 6} < 100, st{1} = sprintf('%s, %.0f%% 2+', st{1}, sets{j, 6}); end
        if o.LayoutCF == 0, st{end+1} = sprintf('CF%d', sets{j, 7}); end
        subtitle(ax, st, 'FontSize', s.FontSizeAnnot - 1);
    end
    cb = colorbar(ax); cb.Layout.Tile = 'east'; cb.Label.String = 'Camera height [m]';
    xlabel(tl, 'x [m]'); ylabel(tl, 'y [m]');
    ttl = sprintf('%s (grey: floor band)', selTitle);
    if any(cellfun(@(c) isfield(c.Mount, 'Cabinet'), cfgs)), ttl = sprintf('%s (brown: cabinets)', selTitle); end
    title(tl, ttl, 'FontSize', s.FontSizeTitle);
    applyThesisStyle(fig);
    exportThesisFigure(fig, fullfile(o.OutDir, [o.Tag '_layouts_' sel]));

    %% Figure 3: per-point occlusion (mean over height layers)
    occ = cell(nP, 1);
    for j = 1:nP
        sp = sets{j, 3};
        [cams, cc] = setupCameras(sets{j, 2}, sp.Cams, sp.Resolution, sp.Focal, sp.FocalWide, ...
            sp.PrincipalPoint, sp.PixelSize);
        [~, occ{j}] = dynamicOcclusion(sp, cams, cc);
    end
    cmax = prctile(vertcat(occ{:}), 98);
    fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 s.FigWidthFull rowH*nRow]);
    tl = tiledlayout(fig, nRow, nCol, 'TileSpacing', 'compact', 'Padding', 'compact');
    for j = 1:nP
        ax = nexttile(tl, tile(j));
        Tg = sets{j, 3}.Target;  ux = unique(Tg(:, 1));  uy = unique(Tg(:, 2));
        M = nan(numel(uy), numel(ux));
        for a = 1:numel(ux)
            for b = 1:numel(uy)
                k = abs(Tg(:, 1) - ux(a)) < 1e-9 & abs(Tg(:, 2) - uy(b)) < 1e-9;
                M(b, a) = mean(occ{j}(k));
            end
        end
        imagesc(ax, ux, uy, M); set(ax, 'YDir', 'normal'); axis(ax, 'equal', 'tight');
        clim(ax, [0 cmax]); colormap(ax, costColormap());
        title(ax, sets{j, 1}, 'FontSize', s.FontSizeAnnot);
        subtitle(ax, sprintf('mean %.0f deg', mean(occ{j})), 'FontSize', s.FontSizeAnnot - 1);
    end
    cb = colorbar(ax); cb.Layout.Tile = 'east'; cb.Label.String = 'Occlusion angle [deg]';
    xlabel(tl, 'x [m]'); ylabel(tl, 'y [m]');
    title(tl, sprintf('Occlusion, %s (mean over heights)', lower(selTitle)), ...
        'FontSize', s.FontSizeTitle);
    applyThesisStyle(fig);
    exportThesisFigure(fig, fullfile(o.OutDir, [o.Tag '_occlusion_' sel]));
end
