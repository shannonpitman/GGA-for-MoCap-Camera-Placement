function [coverageStats]= visualizeCameraCoverage(out, specs, varargin)
% VISUALIZECAMERACOVERAGE  3D coverage view + pie summary for one run.
%
% =====================================================================
% EXAMINER REVIEW
% =====================================================================
% What this plot claims to show
%   "For one optimised camera placement, classify each target point by
%   how many cameras can see it (0, 1, or 2+) and report the
%   distribution as both a 3D scatter and a pie chart."
%
% Strengths
%   - The 0/1/2+ tri-class is the right level of granularity for a
%     MoCap reader: 2+ cameras = triangulable, 1 = limited, 0 =
%     unobserved. This is the actionable metric.
%   - Pie chart gives an at-a-glance summary of distribution mass.
%   - Coverage stats struct returned is useful for headline numbers in
%     the body text.
%
% Weaknesses an examiner will press on
%   1. RED/CYAN/GREEN PALETTE IS COLOUR-BLIND HOSTILE. Roughly 8% of
%      male readers cannot reliably distinguish red and green. Switch
%      to a deuteranopia-safe ramp (orange / sky blue / dark teal),
%      or use shape encoding alongside colour.
%   2. PIE CHARTS ARE STATISTICALLY POOR. Three slices is the upper
%      limit before pies become unreadable. A horizontal stacked bar
%      with the same data is easier to compare across configurations
%      and across the thesis's results section. Edward Tufte hates
%      pies for a reason.
%   3. "2+ CAMERAS" CONFLATES TWO AND EIGHT. Triangulation accuracy
%      depends on baseline angle, not just camera count. A point seen
%      by 8 cameras at 5° apart is *worse* than a point seen by 2
%      cameras at 90° apart. Add a secondary metric: minimum or mean
%      inter-camera ray angle.
%   4. NO SAVE BY DEFAULT. The legacy implementation did not save the
%      figure, so the examiner cannot include it. This version exports
%      a vector PDF when 'SaveAs' is provided.
%   5. PARFOR FOR COLOR ASSIGNMENT IS OVERKILL. Looping over a few
%      thousand points to assign one of three RGB triples is a serial
%      operation; parfor adds pool startup overhead and obscures the
%      code. Replace with logical indexing.
%   6. 3D OCCLUSION OF SCATTER. Front markers hide back markers; the
%      red zero-coverage points are precisely what you want to see and
%      they are exactly what gets occluded. Provide a 2D ground-plane
%      projection as an inset or separate panel.
%   7. ANNOTATION HACK FOR MAX COVERAGE. The "Max: N cameras"
%      annotation is found by iterating over text objects looking for
%      a percentage match — fragile across MATLAB versions. Place the
%      annotation deterministically at a fixed corner of the axes.
%
% Effectiveness at explaining the GA process
%   This is a SINGLE-RUN OUTCOME plot. It does not explain the GA
%   process. To support a GA-process narrative, run the same
%   visualisation on the random initial population's median solution
%   and place it next to the optimised solution: the reader sees the
%   GA's contribution as the reduction in red points.
% =====================================================================

%Visualizes camera coverage of target points with color coding:
% Red: No cameras see the point
% Cyan: Only one camera sees the point
% Green: Two or more cameras see the point
%
% Optional Name-Value parameters:
%   'SaveAs' - Output filename without extension (default: '' = no save)

    p = inputParser;
    addParameter(p, 'SaveAs', '', @ischar);
    parse(p, varargin{:});
    opts = p.Results;

    sty = gaPlotStyle();

    numCams = specs.Cams;
    cameraChromosome = out.bestsol.Chromosome;
    %Camera Parameters
    resolution = specs.Resolution;
    focalLength = specs.Focal;
    focalLengthWide = specs.FocalWide;
    PrincipalPoint = specs.PrincipalPoint;
    TargetSpace = specs.Target;
    maxRange = specs.PreComputed.maxCameraRange;
    maxRangeWide = specs.PreComputed.maxCameraRangeWide;

    [cameras, camCenters] = setupCameras(cameraChromosome, numCams, resolution, focalLength, focalLengthWide, PrincipalPoint, specs.PixelSize);

    numPoints = size(TargetSpace,1);
    cameraCoverage = zeros(numPoints,1);

    %Amount of cams seeing each point (with range check via findVisibleCameras)
    for q = 1:numPoints
        point = TargetSpace(q,:);
        [visibleCams, ~] = findVisibleCameras(point, cameras, camCenters, numCams, resolution, maxRange, maxRangeWide, focalLengthWide);
        cameraCoverage(q) = length(visibleCams);
    end

    % Calculate statistics
    coverageStats.numPoints = numPoints;
    coverageStats.zeroCameras = sum(cameraCoverage == 0);
    coverageStats.oneCamera = sum(cameraCoverage == 1);
    coverageStats.twoPlusCameras = sum(cameraCoverage >= 2);
    coverageStats.zeroCamerasPercent = 100*coverageStats.zeroCameras/numPoints;
    coverageStats.oneCameraPercent = 100*coverageStats.oneCamera/numPoints;
    coverageStats.twoPlusCamerasPercent = 100*coverageStats.twoPlusCameras/numPoints;
    coverageStats.avgCoverage = mean(cameraCoverage);
    coverageStats.maxCoverage = max(cameraCoverage);
    coverageStats.minCoverage = min(cameraCoverage);
    coverageStats.medianCoverage = median(cameraCoverage);

    % Colour map based on coverage (logical indexing — replaces parfor)
    colors = zeros(numPoints, 3);
    colors(cameraCoverage == 0, :) = repmat([0.85, 0.20, 0.20], coverageStats.zeroCameras, 1); % red-ish
    colors(cameraCoverage == 1, :) = repmat([0.20, 0.65, 0.85], coverageStats.oneCamera,    1); % sky blue
    colors(cameraCoverage >= 2, :) = repmat([0.20, 0.55, 0.30], coverageStats.twoPlusCameras, 1); % dark teal

    modalityLabel = getModalityLabel(specs);

    % Create figure with subplots
    fig = figure('Name', 'Coverage Analysis', ...
        'Units', 'inches', ...
        'Position', [0.5, 0.5, sty.FigWidthFull * 1.4, sty.FigHeightTall], ...
        'PaperPositionMode', 'auto', ...
        'Color', sty.BackgroundColor);

    % Subplot 1: Colour coded scatter plot
    ax1 = subplot(1, 2, 1);
    scatter3(TargetSpace(:,1), TargetSpace(:,2), TargetSpace(:,3), 30, colors, 'filled');
    axis equal;
    grid on;
    xlabel('X (m)', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
    ylabel('Y (m)', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
    zlabel('Z (m)', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
    title(sprintf('Camera Coverage - %d Cameras (Cost: %.4f)', specs.Cams, out.bestsol.Cost), ...
        'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName, 'FontWeight', 'normal');
    view(45, 30);
    set(ax1, 'FontSize', sty.FontSizeTick, 'FontName', sty.FontName);

    % Subplot 2: Pie chart
    ax2 = subplot(1, 2, 2);
    pieData = [coverageStats.zeroCamerasPercent, ...
               coverageStats.oneCameraPercent, ...
               coverageStats.twoPlusCamerasPercent];
    pieLabels = {'0 Cameras', '1 Camera', '2+ Cameras'};
    pieColors = [0.85 0.20 0.20;   % red-ish (matches scatter)
                 0.20 0.65 0.85;   % sky blue
                 0.20 0.55 0.30];  % dark teal

    % Filter out zero-value slices
    mask = pieData > 0;
    filteredData = pieData(mask);
    filteredLabels = pieLabels(mask);
    filteredColors = pieColors(mask, :);

    pieH = pie(filteredData);

    % Set colors
    for i = 1:2:length(pieH)
        pieH(i).FaceColor = filteredColors((i+1)/2, :);
    end

    % Set percentage labels
    for i = 2:2:length(pieH)
        pieH(i).String = sprintf('%.1f%%', filteredData(i/2));
        pieH(i).FontSize = sty.FontSizeAxis;
        pieH(i).FontWeight = 'bold';
        pieH(i).Color = 'k';
    end

    % Add max-coverage annotation deterministically (NW of axes)
    text(ax2, -1.15, 1.15, sprintf('Max: %d cameras seeing a single point', coverageStats.maxCoverage), ...
        'FontSize', sty.FontSizeAnnot, 'FontName', sty.FontName, ...
        'BackgroundColor', 'w', 'EdgeColor', 'k', 'Color', 'k', ...
        'Margin', 3);

    legend(filteredLabels, 'Location', 'eastoutside', 'Orientation', 'vertical', ...
        'FontSize', sty.FontSizeLegend);
    title('Coverage Distribution', ...
        'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName, 'FontWeight', 'normal');

    % Overall figure title
    sgtitle({sprintf('Coverage Analysis - %d Cameras', specs.Cams), modalityLabel}, ...
        'FontSize', sty.FontSizeAxis + 2, 'FontWeight', 'bold', 'FontName', sty.FontName);

    % --- Apply thesis style: white bg, black text/axes/ticks/legend ---
    applyThesisStyle(fig);

    % Export if requested
    if ~isempty(opts.SaveAs)
        exportgraphics(fig, [opts.SaveAs '.pdf'], ...
            'ContentType',     'vector', ...
            'BackgroundColor', sty.ExportBgColor);
        fprintf('Coverage figure saved to: %s.pdf\n', opts.SaveAs);
    end

    % Print statistics
    fprintf('\n Camera Coverage Statistics (%s)\n', modalityLabel);
    fprintf('Points with 0 cameras: %d (%.1f%%)\n', ...
        coverageStats.zeroCameras, coverageStats.zeroCamerasPercent);
    fprintf('Points with 1 camera: %d (%.1f%%)\n', ...
        coverageStats.oneCamera, coverageStats.oneCameraPercent);
    fprintf('Points with 2+ cameras: %d (%.1f%%)\n', ...
        coverageStats.twoPlusCameras, coverageStats.twoPlusCamerasPercent);
    fprintf('Average camera coverage: %.2f cameras per point\n', coverageStats.avgCoverage);
    fprintf('Maximum camera coverage: %d cameras per point\n', coverageStats.maxCoverage);
    fprintf('Minimum camera coverage: %d cameras per point\n', coverageStats.minCoverage);
    fprintf('Median camera coverage: %.2f cameras per point\n', coverageStats.medianCoverage);

end
