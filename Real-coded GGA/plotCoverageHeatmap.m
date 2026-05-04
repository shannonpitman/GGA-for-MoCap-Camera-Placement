function plotCoverageHeatmap(varargin)
% PLOTCOVERAGEHEATMAP  Side-by-side 3D coverage heat maps for Uniform vs Normal grid.
%
% Loads the best-performing UAV run for each grid mode from GGA_RunsLog.mat,
% reconstructs camera coverage per target point, and plots 3D scatter heat
% maps coloured by number of visible cameras (blue=0 cameras, red=all cameras).
%
% USAGE:
%   plotCoverageHeatmap()                        % Best UAV run, any camera count
%   plotCoverageHeatmap('NumCameras', 7)         % Best 7-camera UAV run
%   plotCoverageHeatmap('CostFunction', 3)       % Best combined-cost UAV run
%   plotCoverageHeatmap('LogFile', 'mylog.mat')  % Custom log file
%   plotCoverageHeatmap('MarkerSize', 50)        % Larger scatter markers
%
% The function finds the lowest-cost UAV run for each grid mode (Uniform
% and Normal), loads the full result .mat file, and computes per-point
% camera visibility using the same projection logic as visualizeCameraCoverage.
%
% See also: plotGARuns, visualizeCameraCoverage

    %% Parse inputs
    p = inputParser;
    addParameter(p, 'LogFile',      'GGA_RunsLog.mat', @ischar);
    addParameter(p, 'NumCameras',   [],                @isnumeric);
    addParameter(p, 'CostFunction', [],                @isnumeric);
    addParameter(p, 'MarkerSize',   40,                @isnumeric);
    addParameter(p, 'ViewAngle',    [45 30],           @isnumeric);
    parse(p, varargin{:});
    
    opts = p.Results;
    
    %% Load master log
    if ~isfile(opts.LogFile)
        error('Log file not found: %s', opts.LogFile);
    end
    load(opts.LogFile, 'runLog');
    
    % Ensure new fields exist
    requiredFields = {'TargetType', 'GridMode', 'Spacing', 'NumTargetPoints'};
    for f = 1:length(requiredFields)
        if ~isfield(runLog, requiredFields{f})
            [runLog.(requiredFields{f})] = deal(NaN);
        end
    end
    
    %% Find best UAV run for each grid mode
    gridModeNames = {'Uniform Grid', 'Normal Grid'};
    results = cell(1, 2); % {uniform, normal}
    
    for gm = 1:2
        % Filter: UAV (TargetType=1), this grid mode, valid entries
        mask = ([runLog.TargetType] == 1) & ([runLog.GridMode] == gm);
        
        if ~isempty(opts.NumCameras)
            mask = mask & ([runLog.NumCameras] == opts.NumCameras);
        end
        if ~isempty(opts.CostFunction)
            mask = mask & ([runLog.CostFunctionType] == opts.CostFunction);
        end
        
        candidates = runLog(mask);
        
        if isempty(candidates)
            fprintf('No UAV runs found for %s.\n', gridModeNames{gm});
            continue;
        end
        
        % Find best cost
        [bestCost, bestIdx] = min([candidates.BestCost]);
        bestRun = candidates(bestIdx);
        
        % Load the full .mat file
        matFile = bestRun.RunFilename;
        if ~isfile(matFile)
            fprintf('Result file not found: %s (skipping %s)\n', matFile, gridModeNames{gm});
            continue;
        end
        
        loaded = load(matFile, 'saveData');
        results{gm} = loaded.saveData;
        
        fprintf('Loaded %s: %s (Cost=%.4f, %d cameras)\n', ...
            gridModeNames{gm}, matFile, bestCost, bestRun.NumCameras);
    end
    
    %% Check we have at least one result
    hasResults = ~cellfun(@isempty, results);
    if ~any(hasResults)
        error('No matching UAV runs found. Run batch tests first.');
    end
    
    %% Compute coverage and plot
    nPlots = sum(hasResults);
    figure('Name', 'Coverage Heat Map Comparison', 'Position', [50, 50, 700*nPlots, 600]);
    
    plotIdx = 0;
    globalMaxCams = 0;
    
    % First pass: find global max camera count for consistent colormap
    for gm = 1:2
        if hasResults(gm)
            globalMaxCams = max(globalMaxCams, results{gm}.Specifications.Cams);
        end
    end
    
    for gm = 1:2
        if ~hasResults(gm)
            continue;
        end
        
        plotIdx = plotIdx + 1;
        sd = results{gm};
        specs = sd.Specifications;
        chromosome = sd.BestSolution.Chromosome;
        
        % Reconstruct cameras
        numCams = specs.Cams;
        resolution = specs.Resolution;
        focalLength = specs.Focal;
        focalLengthWide = specs.FocalWide;
        PrincipalPoint = specs.PrincipalPoint;
        
        [cameras, ~] = setupCameras(chromosome, numCams, resolution, ...
            focalLength, focalLengthWide, PrincipalPoint);
        
        TargetSpace = specs.Target;
        numPoints = size(TargetSpace, 1);
        
        % Count cameras seeing each point
        cameraCoverage = zeros(numPoints, 1);
        for pt = 1:numPoints
            point = TargetSpace(pt, :);
            visCount = 0;
            for c = 1:numCams
                uv = cameras{c}.project(point);
                if (uv(1) >= 1 && uv(1) <= resolution(1) && ...
                    uv(2) >= 1 && uv(2) <= resolution(2))
                    visCount = visCount + 1;
                end
            end
            cameraCoverage(pt) = visCount;
        end
        
        % Plot 3D scatter with continuous colourmap
        subplot(1, nPlots, plotIdx);
        
        scatter3(TargetSpace(:,1), TargetSpace(:,2), TargetSpace(:,3), ...
            opts.MarkerSize, cameraCoverage, 'filled', 'MarkerFaceAlpha', 0.8);
        
        % Colourmap: blue (0 cameras) -> yellow (mid) -> red (all cameras)
        cmap = buildHeatmapColormap(globalMaxCams);
        colormap(gca, cmap);
        clim([0 globalMaxCams]);
        cb = colorbar;
        cb.Label.String = 'Number of Visible Cameras';
        cb.Label.FontSize = 11;
        
        % Set integer ticks on colourbar
        cb.Ticks = 0:globalMaxCams;
        
        axis equal;
        grid on;
        xlabel('X (m)');
        ylabel('Y (m)');
        zlabel('Z (m)');
        view(opts.ViewAngle);
        
        % Coverage statistics for title
        zeroPct = 100 * sum(cameraCoverage == 0) / numPoints;
        twoPlusPct = 100 * sum(cameraCoverage >= 2) / numPoints;
        avgCov = mean(cameraCoverage);
        
        title({sprintf('%s — %d Cameras (Cost: %.4f)', ...
            gridModeNames{gm}, numCams, sd.BestCost), ...
            sprintf('Avg: %.1f cams/pt | 0-cam: %.1f%% | 2+cam: %.1f%%', ...
            avgCov, zeroPct, twoPlusPct)}, 'FontSize', 11);
    end
    
    % Overall figure title
    if nPlots == 2
        sgtitle('UAV Coverage Heat Map — Uniform vs Normal Grid', ...
            'FontSize', 14, 'FontWeight', 'bold');
    else
        sgtitle('UAV Coverage Heat Map', 'FontSize', 14, 'FontWeight', 'bold');
    end
    
    % Save
    figFilename = 'CoverageHeatmap_UniformVsNormal.png';
    saveas(gcf, figFilename);
    fprintf('Heat map saved to: %s\n', figFilename);
end


function cmap = buildHeatmapColormap(maxCams)
% Build a blue -> cyan -> yellow -> red colourmap with N+1 discrete levels
% so each integer camera count gets a distinct colour.

    nLevels = max(maxCams + 1, 2);
    
    % Define anchor colours
    % 0 cameras = dark blue, mid = yellow/green, max = dark red
    anchors = [
        0.1  0.1  0.6;   % dark blue  (0 cameras)
        0.2  0.5  0.9;   % blue       
        0.2  0.8  0.8;   % cyan       
        0.5  0.9  0.3;   % green-yellow
        0.95 0.9  0.2;   % yellow     
        0.95 0.5  0.1;   % orange     
        0.8  0.1  0.1;   % dark red   (all cameras)
    ];
    
    % Interpolate to nLevels
    nAnchors = size(anchors, 1);
    xi = linspace(1, nAnchors, nLevels);
    cmap = interp1(1:nAnchors, anchors, xi, 'pchip');
    
    % Clamp to [0, 1]
    cmap = max(0, min(1, cmap));
end