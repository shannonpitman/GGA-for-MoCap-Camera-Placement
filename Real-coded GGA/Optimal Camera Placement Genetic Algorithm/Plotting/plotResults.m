function plotResults(out, specs, params, elapsedTime)
% PLOTRESULTS  Single-run convergence + camera-placement visualisation.
%
% =====================================================================
% EXAMINER REVIEW
% =====================================================================
% What this plot claims to show
%   "For one representative GA run: (a) the best/average/top-10 cost
%   converges towards a low asymptote, and (b) the resulting camera
%   placement in 3D world space."
%
% Strengths
%   - Convergence + placement on the same figure tells one complete
%     story about a single run, useful for the methods section.
%   - Log-scale convergence handles wide cost dynamic range.
%   - Top-10 vs population-average vs best curves give a textbook view
%     of GA convergence behaviour.
%
% Weaknesses an examiner will press on
%   1. THIS IS ANECDOTAL EVIDENCE. A single run is not a result; it is
%      an *illustration*. Caption it explicitly as "Example run" so the
%      reader does not mistake it for an aggregate finding. The
%      statistical results live in plotGA_Convergence and the box
%      plots.
%   2. CAMERA PLOT WITHOUT TARGETS. The 3D placement panel shows the
%      cameras' frustums but not the workspace they are observing.
%      The reader cannot judge whether the placement is sensible
%      without seeing what is being captured. Overlay the target
%      space (or at least a wireframe of the workspace) on the same
%      axes.
%   3. TOOLBOX DEPENDENCY. eul2rotm and CentralCamera couple this to
%      the Robotics System Toolbox / Peter Corke RVC toolbox. State
%      this dependency at the top of the file so a reader replicating
%      the work knows what to install.
%   4. PNG SAVES (NOW PDF). Original implementation saved PNG; this
%      version uses exportgraphics + PDF for vector embedding in
%      LaTeX.
%   5. NO TIMESTAMPED SUFFIX BY DEFAULT. The auto-filename includes a
%      datetime, which is good for traceability but bad for build
%      reproducibility (LaTeX includegraphics keys break on re-run).
%      Consider an optional 'Tag' parameter so callers can pin
%      filenames.
%   6. ANNOTATION TEXT ON THE LEFT PANEL OVERLAPS DATA. The
%      `text(...)` call placing "Final Cost / Time" at 0.6·MaxIt may
%      land on the convergence line for short runs. Use a fixed
%      corner annotation or a textbox tied to NW corner instead.
%
% Effectiveness at explaining the GA process
%   The CONVERGENCE panel is good methodology evidence but only at the
%   single-run level. The CAMERA PANEL is a results visualisation but
%   says nothing about how the GA arrived at this layout. To bridge
%   process → result, consider an optional third panel showing the
%   *initial* (random) population's median camera layout next to the
%   final layout — readers see directly what the GA did.
% =====================================================================

currentDateTime = datetime('now');
dateTimeStr = string(currentDateTime, 'yyyyMMdd_HHmmSS');
modalityLabel = getModalityLabel(specs);
sty = gaPlotStyle();

% =====================================================================
% Convergence Figure (left = best cost, right = population evolution)
% =====================================================================
fig1 = figure('Name', 'GA Convergence Analysis', ...
    'Units', 'inches', ...
    'Position', [1, 1, sty.FigWidthFull, sty.FigHeight], ...
    'PaperPositionMode', 'auto', ...
    'Color', sty.BackgroundColor);

% Subplot 1: Best Cost Convergence
ax1 = subplot(1, 2, 1);
semilogy(out.bestcost, 'LineWidth', sty.LineWidth, 'Color', sty.CostFuncColors(1,:));
hold on;
grid on;
xlabel('Iterations', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
ylabel('Best Cost (logarithmic)', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
title(sprintf('Best Cost Convergence - %d Cameras', specs.Cams), ...
    'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName, 'FontWeight', 'normal');
text(0.05*params.MaxIt, max(out.bestcost)*0.5, ...
    sprintf('Final Cost: %.4f\nTime: %.1fsec', out.bestsol.Cost, elapsedTime), ...
    'FontSize', sty.FontSizeAnnot, 'FontName', sty.FontName, ...
    'BackgroundColor', 'w', 'EdgeColor', 'k', 'Color', 'k', 'Margin', 3);
hold off;

% Subplot 2: Average Cost Evolution
ax2 = subplot(1, 2, 2);
semilogy(out.avgcost,        'LineWidth', sty.LineWidth, ...
    'Color', sty.CostFuncColors(2,:), 'DisplayName', 'Population Average');
hold on;
semilogy(out.topTenAvgCost,  'LineWidth', sty.LineWidth, ...
    'Color', sty.CostFuncColors(3,:), 'DisplayName', 'Top 10 Average');
semilogy(out.bestcost,       'LineWidth', sty.LineWidth - 0.2, ...
    'Color', sty.CostFuncColors(1,:), 'LineStyle', '--', 'DisplayName', 'Best Cost');
grid on;
xlabel('Iterations', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
ylabel('Cost (logarithmic)', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
title('Population Cost Evolution', ...
    'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName, 'FontWeight', 'normal');
legend('Location', 'best', 'FontSize', sty.FontSizeLegend);
hold off;

% Force common font / box / tick sizing
for ax = [ax1, ax2]
    set(ax, 'FontSize', sty.FontSizeTick, 'FontName', sty.FontName, ...
        'Box', 'on', 'TickDir', 'out');
end

if specs.warmStart
    warmStr = 'Warm-Started';
else
    warmStr = 'Random Population';
end
sgtitle({sprintf('GA Convergence Analysis - %d Cameras (%s)', specs.Cams, warmStr), modalityLabel}, ...
    'FontSize', sty.FontSizeAxis + 2, 'FontWeight', 'bold', 'FontName', sty.FontName);

% --- Apply thesis style: white bg, black text/axes/ticks/legend ---
applyThesisStyle(fig1);

% Export as vector PDF (replaces legacy PNG)
plotFilename = sprintf('%dCams_Run_%s_convergence', specs.Cams, dateTimeStr);
exportgraphics(fig1, [plotFilename '.pdf'], ...
    'ContentType',     'vector', ...
    'BackgroundColor', sty.ExportBgColor);
fprintf('Convergence plot saved to: %s.pdf\n', plotFilename);

% =====================================================================
% Camera Visualisation Figure
% =====================================================================
% Build cameras with the SAME focal length and resolution that
% setupCameras / the cost function use, so the plotted frustum reflects
% the FOV that actually drives visibility checks.
[cameras, camCenters] = setupCameras(out.bestsol.Chromosome, specs.Cams, ...
    specs.Resolution, specs.Focal, specs.FocalWide, specs.PrincipalPoint);

% Effective range per lens type (matches findVisibleCameras post-fix)
maxRange     = specs.PreComputed.maxCameraRange;
maxRangeWide = specs.PreComputed.maxCameraRangeWide;
focalWide    = specs.FocalWide;

fig2 = figure('Name', 'Camera Configuration', ...
    'Units', 'inches', ...
    'Position', [1, 1, sty.FigWidthFull, sty.FigHeightTall], ...
    'PaperPositionMode', 'auto', ...
    'Color', sty.BackgroundColor);
hold on;

% Plot cameras with frustums clipped at their actual effective range
for i = 1:specs.Cams
    if cameras{i}.f == focalWide
        rng_i = maxRangeWide;
        col_i = [0.85 0.40 0.10]; % wide lens — orange
    else
        rng_i = maxRange;
        col_i = [0.10 0.45 0.75]; % narrow lens — blue
    end
    plotCameraFOV(cameras{i}, camCenters(:,i), rng_i, ...
        'Color', col_i, 'Label', sprintf('cam%d', i));
end

% Overlay target space so the reader can see what's being captured
TargetSpace = specs.Target;
plot3(TargetSpace(:,1), TargetSpace(:,2), TargetSpace(:,3), '.', ...
    'Color', [0.4 0.4 0.4], 'MarkerSize', 3);

axis('equal')
grid('on')
xlabel('X (m)', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName)
ylabel('Y (m)', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName)
zlabel('Z (m)', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName)
title({sprintf('Optimal Camera Placement - Cost: %.4f', out.bestsol.Cost), modalityLabel}, ...
    'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName, 'FontWeight', 'normal')
view(45, 30)
set(gca, 'FontSize', sty.FontSizeTick, 'FontName', sty.FontName);
hold off

% --- Apply thesis style: white bg, black text/axes/ticks ---
applyThesisStyle(fig2);

cameraPlotFilename = sprintf('%dCams_Run_%s_cameras', specs.Cams, dateTimeStr);
exportgraphics(fig2, [cameraPlotFilename '.pdf'], ...
    'ContentType',     'vector', ...
    'BackgroundColor', sty.ExportBgColor);
fprintf('Camera plot saved to: %s.pdf\n', cameraPlotFilename);
end
