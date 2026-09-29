function cfg = analyseConfiguration(varargin)
% ANALYSECONFIGURATION  Analyse and tweak the best GA camera configuration.
%
% Loads the best CF3 (combined cost) result for a given condition, displays
% camera positions and orientations, evaluates all three cost components
% separately, and provides tools for manual position adjustment with
% real-time cost feedback.
%
% USAGE:
%   analyseConfiguration()                              % Best 7C CF3 UAV Uniform
%   analyseConfiguration('NumCameras', 8)               % Best 8C CF3 UAV Uniform
%   analyseConfiguration('TargetType', 2)               % Best UGV result
%   analyseConfiguration('GridMode', 2)                 % Best Normal grid result
%   analyseConfiguration('RunFile', '7Cams_Run_X.mat')  % Specific result file
%
% Returns a struct of function handles for interactive tweaking. Call it
% without an output to just print and plot the loaded configuration.
%
% After loading, use the returned handle to tweak and re-evaluate:
%
%   cfg = analyseConfiguration('NumCameras', 7);
%
%   % View camera positions
%   cfg.printCameras()
%
%   % Manually adjust camera 3 position (e.g., snap to wall)
%   cfg.moveCamera(3, [4.0, -2.0, 2.5], [])     % new position, keep orientation
%   cfg.moveCamera(3, [], [90, -30, 0])           % keep position, new orientation (degrees)
%   cfg.moveCamera(3, [4.0, -2.0, 2.5], [90, -30, 0])  % both
%
%   % Round all positions to nearest 10cm and orientations to nearest 5 degrees
%   cfg.roundAll(0.10, 5)
%
%   % Fix cameras the GA left hanging upside down. The GA never constrains
%   % roll about the optical axis, so optimal solutions routinely invert
%   % cameras. This adds 180 degrees to the roll gene of those cameras,
%   % which spins the image without moving the optical axis, so coverage
%   % and cost are unchanged.
%   cfg.makeUpright()
%   cfg.makeUpright('level')      % stricter: level every horizon (does move cost)
%
%   % Snap orientations to the nearest mountable angle, keeping cameras upright
%   cfg.snapOrientation(15)
%
%   % Snap all cameras to nearest wall (perimeter)
%   cfg.snapToWalls()
%
%   % Re-evaluate all cost components after tweaking
%   cfg.evaluate()
%
%   % Compare original vs tweaked
%   cfg.compareOriginal()
%
%   % Export tweaked configuration for physical setup
%   cfg.exportSetupSheet('7C_setup_sheet.txt')
%
% See also: plotCoverageHeatmap, plotGARuns, batchRunGA

    % Make sure code subfolders are on the path.
    addProjectPaths();

    %% Parse inputs
    % Default log lives under Results/Logs/ (post-restructure layout).
    defaultLog = fullfile(addProjectPaths(), ...
                          'Results', 'Logs', 'GGA_RunsLog.mat');

    p = inputParser;
    addParameter(p, 'LogFile',      defaultLog, @ischar);
    addParameter(p, 'NumCameras',   7,                 @isnumeric);
    addParameter(p, 'TargetType',   1,                 @isnumeric);
    addParameter(p, 'GridMode',     1,                 @isnumeric);
    addParameter(p, 'RunFile',      '',                @ischar);
    addParameter(p, 'Volume',       [-4 4; -4 4; 0 4], @isnumeric);
    parse(p, varargin{:});
    opts = p.Results;

    %% Load result
    if ~isempty(opts.RunFile)
        % Direct file load — accept basename (resolved against Results/) or full path
        runFile = resolveRunPath(opts.RunFile);
        tmp = load(runFile, 'saveData');
        sd = tmp.saveData;
        fprintf('Loaded: %s\n', runFile);
    else
        % Find best CF3 run from log
        load(opts.LogFile, 'runLog');
        
        % Ensure fields exist
        if ~isfield(runLog, 'TargetType'), [runLog.TargetType] = deal(NaN); end
        if ~isfield(runLog, 'GridMode'), [runLog.GridMode] = deal(NaN); end
        
        mask = ([runLog.NumCameras] == opts.NumCameras) & ...
               ([runLog.CostFunctionType] == 3) & ...
               ([runLog.TargetType] == opts.TargetType) & ...
               ([runLog.GridMode] == opts.GridMode);
        filtered = runLog(mask);
        
        if isempty(filtered)
            error('No matching CF3 runs found for %dC TT%d GM%d.', ...
                opts.NumCameras, opts.TargetType, opts.GridMode);
        end
        
        [~, bestIdx] = min([filtered.BestCost]);
        bestRun = filtered(bestIdx);

        bestPath = resolveRunPath(bestRun.RunFilename, bestRun.NumCameras);
        if ~isfile(bestPath)
            error('Result file not found: %s', bestPath);
        end

        tmp = load(bestPath, 'saveData');
        sd = tmp.saveData;
        fprintf('Loaded best CF3: %s (Cost=%.6f)\n', bestPath, bestRun.BestCost);
    end

    %% Extract data
    specs = sd.Specifications;
    numCams = specs.Cams;
    originalChrom = sd.BestSolution.Chromosome;
    currentChrom = originalChrom; % mutable copy

    %% Build output handle (struct of function handles)
    cfg = struct();
    cfg.specs = specs;
    cfg.volume = opts.Volume;
    cfg.originalChromosome = originalChrom;
    cfg.originalCost = sd.BestCost;

    % Print header
    fprintf('\n');
    fprintf('==========================================================\n');
    fprintf('  CONFIGURATION ANALYSIS — %d Cameras\n', numCams);
    fprintf('==========================================================\n');
    ttNames = {'UAV (Full Volume)', 'UGV (Floor Slab)'};
    gmNames = {'Uniform Grid', 'Normal Grid'};
    if isfield(specs, 'TargetType') && specs.TargetType >= 1 && specs.TargetType <= 2
        fprintf('  Target:     %s\n', ttNames{specs.TargetType});
    end
    if isfield(specs, 'TargetMode') && specs.TargetMode >= 1 && specs.TargetMode <= 2
        fprintf('  Grid:       %s\n', gmNames{specs.TargetMode});
    end
    fprintf('  Points:     %d\n', specs.NumPoints);
    fprintf('  Combined Cost: %.6f\n', sd.BestCost);
    fprintf('==========================================================\n\n');

    % --- Print cameras ---
    printCameraTable(currentChrom, numCams, 'OPTIMAL (GA Output)');

    % --- Evaluate all three cost components ---
    fprintf('--- Cost Breakdown (Optimal) ---\n');
    evaluateAllCosts(currentChrom, specs);

    % --- Plot camera configuration ---
    plotCameraConfig(currentChrom, specs, sd.BestCost, 'GA Optimal');

    %% Assign function handles for interactive use
    cfg.printCameras = @() printCameraTable(currentChrom, numCams, 'CURRENT');

    cfg.moveCamera = @(camIdx, newPos, newOriDeg) moveCameraFn(camIdx, newPos, newOriDeg);

    cfg.roundAll = @(posStep, oriStepDeg) roundAllFn(posStep, oriStepDeg);

    cfg.makeUpright = @(varargin) makeUprightFn(varargin{:});

    cfg.snapOrientation = @(oriStepDeg) snapOrientationFn(oriStepDeg);

    cfg.snapToWalls = @() snapToWallsFn();

    cfg.evaluate = @() evaluateCurrentFn();

    cfg.compareOriginal = @() compareOriginalFn();

    cfg.exportSetupSheet = @(filename) exportSetupSheetFn(filename);

    cfg.getChromosome = @() currentChrom;

    cfg.resetToOriginal = @() resetFn();

    %% Nested functions (share currentChrom via closure)

    function moveCameraFn(camIdx, newPos, newOriDeg)
        if camIdx < 1 || camIdx > numCams
            error('Camera index must be between 1 and %d.', numCams);
        end
        chromStart = (camIdx-1)*6 + 1;

        if ~isempty(newPos)
            currentChrom(chromStart:chromStart+2) = newPos;
        end
        if ~isempty(newOriDeg)
            currentChrom(chromStart+3:chromStart+5) = deg2rad(newOriDeg);
        end

        fprintf('Camera %d updated.\n', camIdx);
        printCameraTable(currentChrom, numCams, 'CURRENT (after edit)');
        evaluateAllCosts(currentChrom, specs);
    end

    function roundAllFn(posStep, oriStepDeg)
        % Snapping lives in snapChromosome so this, snapOrientation and the
        % orientationSensitivity study all quantise identically.
        [currentChrom, rep] = snapChromosome(currentChrom, ...
            'OrientationStepDeg', oriStepDeg, 'PositionStep', posStep, ...
            'NumCams', numCams);
        fprintf('Rounded: positions to %.0fcm, orientations to %g degrees.\n', ...
            posStep*100, oriStepDeg);
        fprintf('  Max camera rotation %.2f deg (optical axis moved %.2f deg), max move %.3f m.\n', ...
            max(rep.GeodesicDeg), max(rep.AxisShiftDeg), max(rep.PosShift));
        printCameraTable(currentChrom, numCams, 'CURRENT (after rounding)');
        evaluateAllCosts(currentChrom, specs);
    end

    function makeUprightFn(mode)
        % Roll inverted cameras upright about their own optical axes.
        if nargin < 1 || isempty(mode), mode = 'flip'; end
        [currentChrom, rep] = uprightCameras(currentChrom, ...
            'Mode', mode, 'NumCams', numCams, 'Verbose', true);
        fprintf('%d of %d cameras rolled upright (mode: %s).\n', ...
            rep.NumChanged, numCams, rep.Mode);
        if strcmp(rep.Mode, 'flip')
            fprintf(['  Flip mode changes gamma by 180 deg only, so the optical axis ' ...
                     'and\n  field of view are untouched — the costs below should be ' ...
                     'unchanged.\n']);
        end
        printCameraTable(currentChrom, numCams, ['CURRENT (upright: ' rep.Mode ')']);
        evaluateAllCosts(currentChrom, specs);
    end

    function snapOrientationFn(oriStepDeg)
        % Orientation-only quantisation, positions left alone.
        [currentChrom, rep] = snapChromosome(currentChrom, ...
            'OrientationStepDeg', oriStepDeg, 'NumCams', numCams, 'Verbose', true);
        fprintf('Orientations snapped to the nearest %g degrees.\n', oriStepDeg);
        printCameraTable(currentChrom, numCams, ...
            sprintf('CURRENT (orientations on %g deg grid)', oriStepDeg));
        evaluateAllCosts(currentChrom, specs);
    end

    function snapToWallsFn()
        vol = opts.Volume;
        xMin = vol(1,1); xMax = vol(1,2);
        yMin = vol(2,1); yMax = vol(2,2);

        for c = 1:numCams
            idx = (c-1)*6 + 1;
            x = currentChrom(idx);
            y = currentChrom(idx+1);

            % Find nearest wall
            dists = [abs(x - xMin), abs(x - xMax), abs(y - yMin), abs(y - yMax)];
            [~, wallIdx] = min(dists);

            switch wallIdx
                case 1, currentChrom(idx) = xMin;   % snap to x-min wall
                case 2, currentChrom(idx) = xMax;    % snap to x-max wall
                case 3, currentChrom(idx+1) = yMin;  % snap to y-min wall
                case 4, currentChrom(idx+1) = yMax;  % snap to y-max wall
            end
        end
        wallNames = {'X-min', 'X-max', 'Y-min', 'Y-max'};
        fprintf('All cameras snapped to nearest wall.\n');
        printCameraTable(currentChrom, numCams, 'CURRENT (wall-snapped)');
        evaluateAllCosts(currentChrom, specs);
    end

    function evaluateCurrentFn()
        fprintf('\n--- Cost Breakdown (Current) ---\n');
        evaluateAllCosts(currentChrom, specs);
        plotCameraConfig(currentChrom, specs, NaN, 'Tweaked');
    end

    function compareOriginalFn()
        fprintf('\n');
        fprintf('==========================================================\n');
        fprintf('  COMPARISON: Original vs Tweaked\n');
        fprintf('==========================================================\n\n');

        fprintf('--- ORIGINAL ---\n');
        [origRes, origOcc, origComb] = evaluateAllCosts(originalChrom, specs);

        fprintf('\n--- TWEAKED ---\n');
        [twkRes, twkOcc, twkComb] = evaluateAllCosts(currentChrom, specs);

        fprintf('\n--- DELTA (Tweaked - Original) ---\n');
        fprintf('  Resolution Uncertainty: %+.6f (%+.1f%%)\n', ...
            twkRes - origRes, 100*(twkRes - origRes)/max(origRes, 1e-10));
        fprintf('  Dynamic Occlusion:      %+.6f (%+.1f%%)\n', ...
            twkOcc - origOcc, 100*(twkOcc - origOcc)/max(origOcc, 1e-10));
        fprintf('  Combined Cost:          %+.6f (%+.1f%%)\n', ...
            twkComb - origComb, 100*(twkComb - origComb)/max(origComb, 1e-10));
        fprintf('==========================================================\n\n');

        % Side-by-side camera plot
        figure('Name', 'Original vs Tweaked', 'Position', [50, 50, 1400, 600]);

        subplot(1,2,1);
        plotCamerasInAxis(originalChrom, specs);
        title(sprintf('Original (Cost: %.4f)', origComb));

        subplot(1,2,2);
        plotCamerasInAxis(currentChrom, specs);
        title(sprintf('Tweaked (Cost: %.4f)', twkComb));

        sgtitle(sprintf('Configuration Comparison — %d Cameras', numCams), ...
            'FontSize', 14, 'FontWeight', 'bold');
    end

    function exportSetupSheetFn(filename)
        fid = fopen(filename, 'w');
        fprintf(fid, 'CAMERA SETUP SHEET\n');
        fprintf(fid, '==================\n');
        fprintf(fid, 'Generated: %s\n', char(datetime('now')));
        fprintf(fid, 'Cameras: %d\n\n', numCams);

        [twkRes, twkOcc, twkComb] = evaluateAllCosts(currentChrom, specs);
        fprintf(fid, 'Predicted Costs:\n');
        fprintf(fid, '  Resolution Uncertainty: %.6f\n', twkRes);
        fprintf(fid, '  Dynamic Occlusion:      %.6f\n', twkOcc);
        fprintf(fid, '  Combined:               %.6f\n\n', twkComb);

        fprintf(fid, '%-5s  %-8s %-8s %-8s  %-8s %-8s %-8s  %-8s %-10s %-10s\n', ...
            'Cam', 'X(m)', 'Y(m)', 'Z(m)', 'a(X)', 'b(Y)', 'g(Z)', 'Tilt', 'Mounting', 'Wall');

        fprintf(fid, '%s\n', repmat('-', 1, 95));

        vol = opts.Volume;
        wallTol = 0.15; % 15cm tolerance for wall assignment
        info = cameraOrientationInfo(currentChrom, numCams);

        for c = 1:numCams
            idx = (c-1)*6 + 1;
            pos = currentChrom(idx:idx+2);
            ori = rad2deg(currentChrom(idx+3:idx+5));

            % Determine which wall (if any)
            wallStr = '';
            if abs(pos(1) - vol(1,1)) < wallTol, wallStr = 'X-min';
            elseif abs(pos(1) - vol(1,2)) < wallTol, wallStr = 'X-max';
            elseif abs(pos(2) - vol(2,1)) < wallTol, wallStr = 'Y-min';
            elseif abs(pos(2) - vol(2,2)) < wallTol, wallStr = 'Y-max';
            else, wallStr = 'Interior';
            end

            if info.Degenerate(c)
                mountStr = 'vertical';
            elseif info.Inverted(c)
                mountStr = 'INVERTED';
            else
                mountStr = 'upright';
            end

            fprintf(fid, '%-5d  %-8.2f %-8.2f %-8.2f  %-8.1f %-8.1f %-8.1f  %-8.1f %-10s %-10s\n', ...
                c, pos(1), pos(2), pos(3), ori(1), ori(2), ori(3), ...
                info.TiltDeg(c), mountStr, wallStr);
        end

        nInv = nnz(info.Inverted & ~info.Degenerate);
        if nInv > 0
            fprintf(fid, '\n*** %d camera(s) marked INVERTED must be mounted upside down. ***\n', nInv);
            fprintf(fid, '    Run cfg.makeUpright() before exporting to remove this: it rolls\n');
            fprintf(fid, '    them 180 deg about their own optical axes, which leaves the field\n');
            fprintf(fid, '    of view, the coverage and the cost unchanged.\n');
        end

        fprintf(fid, '\nNotes:\n');
        fprintf(fid, '- Positions in metres relative to room origin\n');
        fprintf(fid, '- Orientations in degrees: intrinsic XYZ Euler angles a, b, g, as\n');
        fprintf(fid, '  R = Rx(a)*Ry(b)*Rz(g). Because Rz is applied last, g is a roll\n');
        fprintf(fid, '  about the optical axis and does not change where the camera looks.\n');
        fprintf(fid, '- Tilt is the roll measured from level: 0 = horizon level in frame,\n');
        fprintf(fid, '  +/-180 = upside down.\n');
        fprintf(fid, '- Wall assignment based on %.0fcm tolerance\n', wallTol*100);
        fprintf(fid, '- Use OptiTrack Motive reprojection overlay to fine-align\n');

        fclose(fid);
        fprintf('Setup sheet exported: %s\n', filename);
    end

    function resetFn()
        currentChrom = originalChrom;
        fprintf('Reset to original GA solution.\n');
        printCameraTable(currentChrom, numCams, 'RESET TO ORIGINAL');
    end
end


%% ====================================================================
%  LOCAL HELPER FUNCTIONS
%  ====================================================================

function printCameraTable(chrom, numCams, label)
    % Tilt is the roll about the optical axis measured from level: 0 means
    % the horizon is level in frame, +/-180 means the camera is inverted.
    % The GA does not constrain it, so it is worth seeing at a glance.
    info = cameraOrientationInfo(chrom, numCams);

    fprintf('\n  %-5s  %-8s %-8s %-8s  %-8s %-8s %-8s  %-8s %-9s\n', ...
        'Cam', 'X(m)', 'Y(m)', 'Z(m)', 'α(X)°', 'β(Y)°', 'γ(Z)°', 'Tilt°', 'Status');
    fprintf('  %s\n', repmat('-', 1, 82));

    for c = 1:numCams
        idx = (c-1)*6 + 1;
        pos = chrom(idx:idx+2);
        ori = rad2deg(chrom(idx+3:idx+5));
        if info.Degenerate(c)
            status = 'vertical';
        elseif info.Inverted(c)
            status = 'INVERTED';
        else
            status = 'upright';
        end
        fprintf('  %-5d  %-8.2f %-8.2f %-8.2f  %-8.1f %-8.1f %-8.1f  %-8.1f %-9s\n', ...
            c, pos(1), pos(2), pos(3), ori(1), ori(2), ori(3), ...
            info.TiltDeg(c), status);
    end

    nInv = nnz(info.Inverted & ~info.Degenerate);
    if nInv > 0
        fprintf(['  %d camera(s) upside down. cfg.makeUpright() rolls them 180° about\n' ...
                 '  their optical axes, which leaves coverage and cost unchanged.\n'], nInv);
    end
    fprintf('  [%s]\n\n', label);
end


function [resCost, occCost, combCost] = evaluateAllCosts(chrom, specs)
    % Evaluate all three cost components
    resCost  = resUncertaintyCost(chrom, specs);
    occCost  = dynamicOcclusionCost(chrom, specs);
    combCost = combinedCostFunction(chrom, specs);

    fprintf('  Resolution Uncertainty: %.6f\n', resCost);
    fprintf('  Dynamic Occlusion:      %.6f\n', occCost);
    fprintf('  Combined Cost:          %.6f\n', combCost);
end


function plotCameraConfig(chrom, specs, cost, label)
    numCams = specs.Cams;

    % Build cameras with the actual focal length, resolution and principal
    % point that the cost function uses, so the drawn FOV matches the FOV
    % seen by findVisibleCameras / computePointUncertainty.
    [cameras, camCenters] = setupCameras(chrom, numCams, specs.Resolution, ...
        specs.Focal, specs.FocalWide, specs.PrincipalPoint, specs.PixelSize);

    figure('Name', sprintf('Camera Config — %s', label), ...
        'Position', [100, 100, 800, 600]);
    hold on;

    drawCamerasWithFOV(cameras, camCenters, specs);
    if isfield(specs, 'Target')
        T = specs.Target;
        plot3(T(:,1), T(:,2), T(:,3), '.', ...
            'Color', [0.4 0.4 0.4], 'MarkerSize', 3);
    end

    axis equal; grid on;
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    if ~isnan(cost)
        title(sprintf('%s — %d Cameras (Cost: %.4f)', label, numCams, cost));
    else
        title(sprintf('%s — %d Cameras', label, numCams));
    end
    view(45, 30);
    hold off;
end


function plotCamerasInAxis(chrom, specs)
    numCams = specs.Cams;

    [cameras, camCenters] = setupCameras(chrom, numCams, specs.Resolution, ...
        specs.Focal, specs.FocalWide, specs.PrincipalPoint, specs.PixelSize);

    hold on;
    drawCamerasWithFOV(cameras, camCenters, specs);
    if isfield(specs, 'Target')
        T = specs.Target;
        plot3(T(:,1), T(:,2), T(:,3), '.', ...
            'Color', [0.4 0.4 0.4], 'MarkerSize', 3);
    end
    axis equal; grid on;
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    view(45, 30);
    hold off;
end


function drawCamerasWithFOV(cameras, camCenters, specs)
% Draw every camera's true FOV pyramid clipped at its effective range.
    numCams      = specs.Cams;
    maxRange     = specs.PreComputed.maxCameraRange;
    maxRangeWide = specs.PreComputed.maxCameraRangeWide;
    focalWide    = specs.FocalWide;

    for i = 1:numCams
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
end