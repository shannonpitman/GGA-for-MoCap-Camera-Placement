function plotLayoutAndCostField(varargin)
%PLOTLAYOUTANDCOSTFIELD  Optimised vs manual, side by side, three rows.
%
%   Column 1 = optimised GA rig, column 2 = manually posed rig.
%   Row 1 = camera poses, row 2 = resolution-uncertainty field,
%   row 3 = dynamic-occlusion field. Colour scales are shared per row so
%   the two columns are directly comparable.
%
%   USAGE
%       plotLayoutAndCostField('TargetType', 1, 'GridMode', 1)

    p = inputParser;
    addParameter(p, 'TargetType',   1,   @isnumeric);
    addParameter(p, 'GridMode',     1,   @isnumeric);
    addParameter(p, 'Spacing',      1.0, @isnumeric);
    addParameter(p, 'NumCameras',   7,   @isnumeric);
    addParameter(p, 'CostFunction', 3,   @isnumeric);
    addParameter(p, 'ViewAngle',    [45 25], @isnumeric);
    addParameter(p, 'PyramidScale', 1.4, @isnumeric);
    addParameter(p, 'ZStretch',     1.0, @isnumeric);   % >1 exaggerates height
    addParameter(p, 'SlideFonts',   false, @islogical); % large text for projection
    addParameter(p, 'MarkerSize',   16,  @isnumeric);
    addParameter(p, 'FigHeight',    6.6, @isnumeric);
    addParameter(p, 'SaveAs',       '',  @ischar);
    parse(p, varargin{:});
    o = p.Results;

    sty   = gaPlotStyle();
    ttStr = sty.TargetNames{o.TargetType};
    gmStr = sty.GridNames{o.GridMode};

    [gaChrom, specs, gaCost] = localBestGARun(o);
    [optiChrom, optiIDs] = buildOptiTrackChromosome();
    optiSpecs      = specs;
    optiSpecs.Cams = 7;
    oc = evaluateOptiTrackCost('TargetType', o.TargetType, ...
                               'GridMode',   o.GridMode, ...
                               'Spacing',    o.Spacing);
    optiCost = oc.CF3;

    [uncGA,   occGA]   = localPerPointCosts(gaChrom,   specs);
    [uncOpti, occOpti] = localPerPointCosts(optiChrom, optiSpecs);

    fig = figure('Name', 'Layout and cost field', 'Units', 'inches', ...
        'Position', [0.5 0.5 13.0 o.FigHeight], ...
        'PaperPositionMode', 'auto', 'Color', sty.BackgroundColor);
    tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    T      = specs.Target;
    uncMax = max([uncGA; uncOpti]);
    occMax = max([occGA; occOpti]);
    if o.SlideFonts
        names = {sprintf('Optimised GA Rig: J = %.3f', gaCost), ...
                 sprintf('Manually Posed Rig: J = %.2f', optiCost)};
    else
        names = {sprintf('Optimised GA Rig  (J = %.4f)', gaCost), ...
                 sprintf('Manually Posed Rig  (J = %.4f)', optiCost)};
    end
    chroms = {gaChrom, optiChrom};
    camIDs = {[], optiIDs};   % label manual cameras with their Motive IDs
    specsC = {specs, optiSpecs};
    uncs   = {uncGA, uncOpti};
    occs   = {occGA, occOpti};

    axAll = gobjects(1,6);
    for k = 1:2
        base = (k-1)*3;

        ax = nexttile(tl, base+1);  axAll(base+1) = ax;
        [hN, hW] = localDrawCameras(ax, chroms{k}, specsC{k}, o, camIDs{k});
        title(ax, names{k}, 'FontSize', sty.FontSizeAxis, ...
            'FontName', sty.FontName, 'FontWeight', 'bold', 'Color', 'k');

        ax = nexttile(tl, base+2); axAll(base+2) = ax;
        localScatter(ax, T, uncs{k}, [0 uncMax], o);
        if o.SlideFonts
            tStr = sprintf('Mean uncertainty: %.3f m', mean(uncs{k}));
        else
            tStr = sprintf('Resolution uncertainty: mean %.3f m', mean(uncs{k}));
        end
        title(ax, tStr, ...
            'FontSize', sty.FontSizeAnnot, 'FontName', sty.FontName, ...
            'FontWeight', 'normal', 'Color', 'k');

        ax = nexttile(tl, base+3); axAll(base+3) = ax;
        localScatter(ax, T, occs{k}, [0 occMax], o);
        if o.SlideFonts
            tStr = sprintf('Mean occlusion: %.0f%s', mean(occs{k}), char(176));
        else
            tStr = sprintf('Dynamic occlusion: mean %.0f deg', mean(occs{k}));
        end
        title(ax, tStr, ...
            'FontSize', sty.FontSizeAnnot, 'FontName', sty.FontName, ...
            'FontWeight', 'normal', 'Color', 'k');
    end

    if o.SlideFonts
        lensLbl = {sprintf('Narrow (%.1f mm)', specs.Focal*1e3), ...
                   sprintf('Wide (%.1f mm)',   specs.FocalWide*1e3)};
    else
        lensLbl = {sprintf('Narrow lens (%.1f mm)', specs.Focal*1e3), ...
                   sprintf('Wide lens (%.1f mm)',   specs.FocalWide*1e3)};
    end
    lg = legend(axAll(1), [hN hW], lensLbl, ...
        'Orientation', 'horizontal', 'FontSize', sty.FontSizeAnnot, ...
        'FontName', sty.FontName, 'Location', 'northoutside');

    cb = colorbar(axAll(5), 'Location', 'eastoutside');
    cb.Label.String = 'Uncertainty (m)';
    cb = colorbar(axAll(6), 'Location', 'eastoutside');
    cb.Label.String = 'Occlusion (deg)';

    applyThesisStyle(fig);

    % Large type for a projected slide: set last, because both setting an
    % axes FontSize and applyThesisStyle would otherwise rescale labels.
    if o.SlideFonts
        for q = 1:6
            a = axAll(q);
            set(a, 'FontSize', 15, 'XTick', [-4 0 4], 'YTick', [-4 0 4], ...
                   'ZTick', [0 2 4]);
            a.XLabel.FontSize = 16;  a.YLabel.FontSize = 16;  a.ZLabel.FontSize = 16;
            a.Title.FontSize  = 19;
        end
        lg.FontSize = 18;
        for c = findall(fig, 'Type', 'colorbar')'
            c.FontSize = 15;  c.Label.FontSize = 17;
        end
    end

    outName = o.SaveAs;
    if isempty(outName)
        outName = sprintf('LayoutAndCostField_%s_%dC', ttStr, o.NumCameras);
    end
    exportThesisFigure(fig, outName, 'Background', sty.ExportBgColor, 'Quiet', true);
    fprintf('Saved: %s.{pdf,png}\n', outName);
end


function localScatter(ax, T, vals, clim, o)
    scatter3(ax, T(:,1), T(:,2), T(:,3), o.MarkerSize, vals, 'filled');
    set(ax, 'CLim', clim); colormap(ax, costColormap(256));
    daspect(ax, [1 1 1/o.ZStretch]); grid(ax, 'on'); view(ax, [45 25]);
    xlabel(ax, 'X (m)'); ylabel(ax, 'Y (m)'); zlabel(ax, 'Z (m)');
end


function [hN, hW] = localDrawCameras(ax, chrom, specs, o, ids)
    numCams = specs.Cams;
    [cameras, C] = setupCameras(chrom, numCams, specs.Resolution, ...
        specs.Focal, specs.FocalWide, specs.PrincipalPoint, specs.PixelSize);
    narrowCol = [0.10 0.45 0.75];  wideCol = [0.85 0.40 0.10];
    hold(ax, 'on');
    T  = specs.Target;
    bb = [min(T,[],1); max(T,[],1)];
    v  = [bb(1,1) bb(1,2) bb(1,3); bb(2,1) bb(1,2) bb(1,3); ...
          bb(2,1) bb(2,2) bb(1,3); bb(1,1) bb(2,2) bb(1,3); ...
          bb(1,1) bb(1,2) bb(2,3); bb(2,1) bb(1,2) bb(2,3); ...
          bb(2,1) bb(2,2) bb(2,3); bb(1,1) bb(2,2) bb(2,3)];
    e = [1 2;2 3;3 4;4 1;5 6;6 7;7 8;8 5;1 5;2 6;3 7;4 8];
    for q = 1:size(e,1)
        plot3(ax, v(e(q,:),1), v(e(q,:),2), v(e(q,:),3), '-', ...
            'Color', [0.35 0.35 0.35], 'LineWidth', 0.6);
    end
    plot3(ax, T(:,1), T(:,2), T(:,3), '.', 'Color', [0.6 0.6 0.6], ...
        'MarkerSize', 2);
    pyrLen  = o.PyramidScale;
    pyrHalf = 0.5 * pyrLen / 2;
    for i = 1:numCams
        if cameras{i}.f == specs.FocalWide, col = wideCol; else, col = narrowCol; end
        R = cameras{i}.T.rotm;
        localPosePyramid(ax, C(:,i), R(:,3), R(:,1), R(:,2), pyrLen, pyrHalf, col);
        plot3(ax, C(1,i), C(2,i), C(3,i), 'o', 'MarkerSize', 4, ...
            'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        if ~isempty(ids)
            text(ax, C(1,i), C(2,i), C(3,i) + 0.45, sprintf('%d', ids(i)), ...
                'FontSize', 8 + 6*o.SlideFonts, 'FontWeight', 'bold', 'Color', col, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        end
    end
    % Legend proxies (not drawn)
    hN = patch(ax, NaN, NaN, NaN, narrowCol, 'FaceAlpha', 0.5, 'EdgeColor', narrowCol);
    hW = patch(ax, NaN, NaN, NaN, wideCol,   'FaceAlpha', 0.5, 'EdgeColor', wideCol);
    daspect(ax, [1 1 1/o.ZStretch]); grid(ax, 'on'); view(ax, o.ViewAngle);
    xlabel(ax, 'X (m)'); ylabel(ax, 'Y (m)'); zlabel(ax, 'Z (m)');
    hold(ax, 'off');
end


function [chrom, specs, cost] = localBestGARun(o)
    root = addProjectPaths();
    S = load(fullfile(root, 'Results', 'Logs', 'GGA_RunsLog.mat'), 'runLog');
    r = S.runLog;
    mask = ([r.NumCameras] == o.NumCameras) & ([r.CostFunctionType] == o.CostFunction) & ...
           ([r.TargetType] == o.TargetType) & ([r.GridMode] == o.GridMode) & ...
           (abs([r.Spacing] - o.Spacing) < 1e-6);
    cand = r(mask);
    [cost, idx] = min([cand.BestCost]);
    L = load(resolveRunPath(cand(idx).RunFilename, cand(idx).NumCameras), 'saveData');
    chrom = L.saveData.BestSolution.Chromosome;
    specs = L.saveData.Specifications;
end


function [unc, occ] = localPerPointCosts(chrom, specs)
% Per-point values straight from the cost functions, so the plotted field
% uses exactly the same visibility and penalties as the GA.
    [cameras, camCenters] = setupCameras(chrom, specs.Cams, specs.Resolution, ...
        specs.Focal, specs.FocalWide, specs.PrincipalPoint, specs.PixelSize);
    [~, unc] = resUncertainty(specs, cameras, camCenters);
    [~, occ] = dynamicOcclusion(specs, cameras, camCenters);
end


function localPosePyramid(ax, apex, axisDir, rightV, upV, len, halfBase, col)
% Four-sided pyramid: apex at the camera centre, base at apex + axis*len.
    apex    = apex(:);
    axisDir = axisDir(:) / norm(axisDir);
    rightV  = rightV(:)  / norm(rightV);
    upV     = upV(:)     / norm(upV);
    c = apex + axisDir * len;
    b = [c + halfBase*( rightV + upV), ...
         c + halfBase*(-rightV + upV), ...
         c + halfBase*(-rightV - upV), ...
         c + halfBase*( rightV - upV)];
    for q = 1:4
        r = mod(q, 4) + 1;
        patch(ax, 'XData', [apex(1) b(1,q) b(1,r)], ...
                  'YData', [apex(2) b(2,q) b(2,r)], ...
                  'ZData', [apex(3) b(3,q) b(3,r)], ...
                  'FaceColor', col, 'FaceAlpha', 0.28, ...
                  'EdgeColor', col, 'LineWidth', 0.5);
    end
    patch(ax, 'XData', b(1,:), 'YData', b(2,:), 'ZData', b(3,:), ...
          'FaceColor', col, 'FaceAlpha', 0.20, ...
          'EdgeColor', col, 'LineWidth', 0.5);
end
