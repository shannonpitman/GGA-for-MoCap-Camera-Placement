classdef ConfigAnalyser < handle
% CONFIGANALYSER  Interactive camera configuration analysis and tweaking.
%
% Loads the best CF3 result, evaluates all three cost components, provides
% wall-snapping, rounding, per-camera editing, and FOV visualisation.
%
% USAGE:
%   cfg = ConfigAnalyser(7, 1, 1);          % 7 cameras, UAV, Uniform grid
%   cfg = ConfigAnalyser(7, 2, 1);          % 7 cameras, UGV, Uniform grid
%   cfg = ConfigAnalyser('File', '7Cams_Run_XXX.mat');  % Specific file
%
%   cfg.printCameras()                      % Print position/orientation table
%   cfg.evaluate()                          % Evaluate all 3 cost components
%   cfg.snapToWalls()                       % Snap each camera to nearest wall
%   cfg.roundAll(0.10, 5)                   % Round pos to 10cm, ori to 5 deg
%   cfg.moveCamera(3, [4 -2 2.5], [90 -30 0])  % Move camera 3
%   cfg.compareOriginal()                   % Side-by-side cost comparison
%   cfg.plotWithFOV()                       % 3D plot with FOV frustums
%   cfg.plotWithFOV('FrustumLength', 3)     % Custom frustum depth
%   cfg.exportSetupSheet('setup.txt')       % Export for lab use
%   cfg.reset()                             % Revert to GA original

    properties
        specs           % Hardware specs struct
        volume          % Capture/workspace volume [3x2]
        wallBounds      % Room walls where cameras can mount [3x2]
        originalChrom   % Original GA chromosome
        currentChrom    % Mutable working chromosome
        originalCost    % Original combined cost
        numCams         % Number of cameras
        sourceFile      % Source .mat filename
    end

    methods
        function obj = ConfigAnalyser(varargin)
            % Make sure code subfolders are on the path.
            addProjectPaths();

            p = inputParser;
            addOptional(p, 'NumCameras', 7, @isnumeric);
            addOptional(p, 'TargetType', 1, @isnumeric);
            addOptional(p, 'GridMode', 1, @isnumeric);
            addParameter(p, 'File', '', @ischar);
            % Default log lives under Results/Logs/
            defaultLog = fullfile(fileparts(mfilename('fullpath')), ...
                                  'Results', 'Logs', 'GGA_RunsLog.mat');
            addParameter(p, 'LogFile', defaultLog, @ischar);
            addParameter(p, 'Volume', [-4 4; -4 4; 0 4], @isnumeric);
            addParameter(p, 'WallBounds', [-5 5; -4.5 4.5; 0 4.8], @isnumeric);
            parse(p, varargin{:});

            obj.volume = p.Results.Volume;
            obj.wallBounds = p.Results.WallBounds;

            if ~isempty(p.Results.File)
                % Load specific file (basename or full path both accepted)
                runFile = resolveRunPath(p.Results.File);
                tmp = load(runFile, 'saveData');
                sd = tmp.saveData;
                obj.sourceFile = runFile;
            else
                % Find best CF3 from log
                load(p.Results.LogFile, 'runLog');
                if ~isfield(runLog, 'TargetType'), [runLog.TargetType] = deal(NaN); end
                if ~isfield(runLog, 'GridMode'), [runLog.GridMode] = deal(NaN); end

                mask = ([runLog.NumCameras] == p.Results.NumCameras) & ...
                       ([runLog.CostFunctionType] == 3) & ...
                       ([runLog.TargetType] == p.Results.TargetType) & ...
                       ([runLog.GridMode] == p.Results.GridMode);
                filtered = runLog(mask);

                if isempty(filtered)
                    error('No matching CF3 runs found for %dC TT%d GM%d.', ...
                        p.Results.NumCameras, p.Results.TargetType, p.Results.GridMode);
                end

                [~, bestIdx] = min([filtered.BestCost]);
                bestRun = filtered(bestIdx);

                bestPath = resolveRunPath(bestRun.RunFilename, bestRun.NumCameras);
                if ~isfile(bestPath)
                    error('Result file not found: %s', bestPath);
                end

                tmp = load(bestPath, 'saveData');
                sd = tmp.saveData;
                obj.sourceFile = bestPath;
            end

            obj.specs = sd.Specifications;
            obj.numCams = obj.specs.Cams;
            obj.originalChrom = sd.BestSolution.Chromosome;
            obj.currentChrom = obj.originalChrom;
            obj.originalCost = sd.BestCost;

            % Print summary
            ttNames = {'UAV (Full Volume)', 'UGV (Floor Slab)'};
            gmNames = {'Uniform Grid', 'Normal Grid'};
            fprintf('\n=== ConfigAnalyser: %d Cameras ===\n', obj.numCams);
            fprintf('Source: %s\n', obj.sourceFile);
            if isfield(obj.specs, 'TargetType') && obj.specs.TargetType >= 1 && obj.specs.TargetType <= 2
                fprintf('Target: %s\n', ttNames{obj.specs.TargetType});
            end
            if isfield(obj.specs, 'TargetMode') && obj.specs.TargetMode >= 1 && obj.specs.TargetMode <= 2
                fprintf('Grid:   %s\n', gmNames{obj.specs.TargetMode});
            end
            fprintf('Points: %d\n', obj.specs.NumPoints);
            fprintf('Original Combined Cost: %.6f\n', obj.originalCost);
            fprintf('==================================\n\n');

            obj.printCameras();
            obj.evaluate();
        end

        %% ============================================================
        %  PRINT CAMERA TABLE
        %  ============================================================
        function printCameras(obj)
            fprintf('\n  %-5s %-6s  %-8s %-8s %-8s  %-8s %-8s %-8s\n', ...
                'Cam', 'Lens', 'X(m)', 'Y(m)', 'Z(m)', 'Roll°', 'Pitch°', 'Yaw°');
            fprintf('  %s\n', repmat('-', 1, 68));

            for c = 1:obj.numCams
                idx = (c-1)*6 + 1;
                pos = obj.currentChrom(idx:idx+2);
                ori = rad2deg(obj.currentChrom(idx+3:idx+5));

                if mod(c,2) == 0
                    lensStr = 'Wide';
                else
                    lensStr = 'Narr';
                end

                fprintf('  %-5d %-6s  %-8.2f %-8.2f %-8.2f  %-8.1f %-8.1f %-8.1f\n', ...
                    c, lensStr, pos(1), pos(2), pos(3), ori(1), ori(2), ori(3));
            end
            fprintf('\n');
        end

        %% ============================================================
        %  EVALUATE ALL COST COMPONENTS
        %  ============================================================
        function [resCost, occCost, combCost] = evaluate(obj)
            resCost  = resUncertaintyCost(obj.currentChrom, obj.specs);
            occCost  = dynamicOcclusionCost(obj.currentChrom, obj.specs);
            combCost = combinedCostFunction(obj.currentChrom, obj.specs);

            fprintf('  Cost Breakdown:\n');
            fprintf('    Resolution Uncertainty: %.6f\n', resCost);
            fprintf('    Dynamic Occlusion:      %.6f\n', occCost);
            fprintf('    Combined Cost:          %.6f\n\n', combCost);
        end

        %% ============================================================
        %  MOVE A SINGLE CAMERA
        %  ============================================================
        function moveCamera(obj, camIdx, newPos, newOriDeg)
        % moveCamera(camIdx, [x y z], [roll pitch yaw])
        % Pass [] for either argument to keep current value.
            if camIdx < 1 || camIdx > obj.numCams
                error('Camera index must be between 1 and %d.', obj.numCams);
            end
            idx = (camIdx-1)*6 + 1;

            if ~isempty(newPos)
                obj.currentChrom(idx:idx+2) = newPos;
            end
            if ~isempty(newOriDeg)
                obj.currentChrom(idx+3:idx+5) = deg2rad(newOriDeg);
            end

            fprintf('Camera %d updated.\n', camIdx);
            obj.printCameras();
            obj.evaluate();
        end

        %% ============================================================
        %  SNAP ALL CAMERAS TO NEAREST WALL
        %  ============================================================
        function snapToWalls(obj)
        % Snaps each camera's x or y coordinate to the nearest room wall.
        % Uses wallBounds (room walls), NOT volume (capture zone).
        % Only the single closest coordinate is moved per camera.
            wb = obj.wallBounds;
            vol = obj.volume;
            xWallMin = wb(1,1); xWallMax = wb(1,2);
            yWallMin = wb(2,1); yWallMax = wb(2,2);
            wallNames = {'X-min', 'X-max', 'Y-min', 'Y-max'};

            fprintf('  Snapping cameras to nearest room wall:\n');
            fprintf('  (Room walls: X=[%.1f, %.1f], Y=[%.1f, %.1f])\n', ...
                xWallMin, xWallMax, yWallMin, yWallMax);
            fprintf('  (Capture volume: X=[%.1f, %.1f], Y=[%.1f, %.1f])\n\n', ...
                vol(1,1), vol(1,2), vol(2,1), vol(2,2));

            for c = 1:obj.numCams
                idx = (c-1)*6 + 1;
                x = obj.currentChrom(idx);
                y = obj.currentChrom(idx+1);

                % Distance to each room wall
                dists = [abs(x - xWallMin), abs(x - xWallMax), ...
                         abs(y - yWallMin), abs(y - yWallMax)];
                [minDist, wallIdx] = min(dists);

                switch wallIdx
                    case 1, obj.currentChrom(idx)   = xWallMin;
                    case 2, obj.currentChrom(idx)   = xWallMax;
                    case 3, obj.currentChrom(idx+1) = yWallMin;
                    case 4, obj.currentChrom(idx+1) = yWallMax;
                end

                fprintf('    Cam %d -> %s wall (moved %.2fm)\n', ...
                    c, wallNames{wallIdx}, minDist);
            end
            fprintf('\n');
            obj.printCameras();
            obj.evaluate();
        end

        %% ============================================================
        %  ROUND ALL POSITIONS AND ORIENTATIONS
        %  ============================================================
        function roundAll(obj, posStep, oriStepDeg)
        % roundAll(positionStep_m, orientationStep_deg)
        %   e.g. roundAll(0.10, 5) rounds to 10cm and 5 degrees
            oriStepRad = deg2rad(oriStepDeg);
            for c = 1:obj.numCams
                idx = (c-1)*6 + 1;
                obj.currentChrom(idx:idx+2) = ...
                    round(obj.currentChrom(idx:idx+2) / posStep) * posStep;
                obj.currentChrom(idx+3:idx+5) = ...
                    round(obj.currentChrom(idx+3:idx+5) / oriStepRad) * oriStepRad;
            end
            fprintf('Rounded: positions to %.0fcm, orientations to %d degrees.\n\n', ...
                posStep*100, oriStepDeg);
            obj.printCameras();
            obj.evaluate();
        end

        %% ============================================================
        %  COMPARE ORIGINAL vs TWEAKED
        %  ============================================================
        function compareOriginal(obj)
            fprintf('\n=== COMPARISON: Original vs Tweaked ===\n\n');

            fprintf('--- ORIGINAL ---\n');
            origRes  = resUncertaintyCost(obj.originalChrom, obj.specs);
            origOcc  = dynamicOcclusionCost(obj.originalChrom, obj.specs);
            origComb = combinedCostFunction(obj.originalChrom, obj.specs);
            fprintf('  Res: %.6f | Occ: %.6f | Comb: %.6f\n\n', origRes, origOcc, origComb);

            fprintf('--- TWEAKED ---\n');
            twkRes  = resUncertaintyCost(obj.currentChrom, obj.specs);
            twkOcc  = dynamicOcclusionCost(obj.currentChrom, obj.specs);
            twkComb = combinedCostFunction(obj.currentChrom, obj.specs);
            fprintf('  Res: %.6f | Occ: %.6f | Comb: %.6f\n\n', twkRes, twkOcc, twkComb);

            fprintf('--- DELTA ---\n');
            fprintf('  Res:  %+.6f (%+.1f%%)\n', twkRes-origRes, 100*(twkRes-origRes)/max(origRes,1e-10));
            fprintf('  Occ:  %+.6f (%+.1f%%)\n', twkOcc-origOcc, 100*(twkOcc-origOcc)/max(origOcc,1e-10));
            fprintf('  Comb: %+.6f (%+.1f%%)\n', twkComb-origComb, 100*(twkComb-origComb)/max(origComb,1e-10));
            fprintf('==========================================\n\n');

            % Side-by-side figure
            figure('Name', 'Original vs Tweaked', 'Position', [50 50 1400 600]);
            subplot(1,2,1);
            obj.plotCamerasInAxis(obj.originalChrom);
            title(sprintf('Original (Cost: %.4f)', origComb));
            subplot(1,2,2);
            obj.plotCamerasInAxis(obj.currentChrom);
            title(sprintf('Tweaked (Cost: %.4f)', twkComb));
            sgtitle(sprintf('Configuration Comparison — %d Cameras', obj.numCams), ...
                'FontSize', 14, 'FontWeight', 'bold');
        end

        %% ============================================================
        %  PLOT WITH FOV FRUSTUMS
        %  ============================================================
        function plotWithFOV(obj, varargin)
        % plotWithFOV()                    % Uses hardware range per lens
        % plotWithFOV('ScaleFOV', 0.5)     % Show half the effective range
        % plotWithFOV('ShowTarget', false) % Hide target points
            ip = inputParser;
            addParameter(ip, 'ScaleFOV', 1.0, @isnumeric);  % fraction of effective range
            addParameter(ip, 'ShowTarget', true, @islogical);
            parse(ip, varargin{:});
            fovScale = ip.Results.ScaleFOV;
            showTarget = ip.Results.ShowTarget;

            figure('Name', 'Camera FOV Visualisation', ...
                'Position', [50 50 1000 800]);
            hold on;

            resolution = obj.specs.Resolution;
            fNarrow = obj.specs.Focal;
            fWide = obj.specs.FocalWide;
            pixSize = obj.specs.PixelSize;
            rangeNarrow = obj.specs.Range;      % effective range [m]
            rangeWide   = obj.specs.RangeWide;  % effective range [m]

            % Colour scheme: narrow=blue, wide=red
            colNarrow = [0.2 0.4 0.9];
            colWide   = [0.9 0.3 0.2];

            fprintf('  FOV Visualisation:\n');
            fprintf('    Narrow lens: f=%.1fmm, range=%.0fm\n', fNarrow*1000, rangeNarrow);
            fprintf('    Wide lens:   f=%.1fmm, range=%.0fm\n', fWide*1000, rangeWide);
            if fovScale ~= 1.0
                fprintf('    Scale: %.0f%% of effective range\n', fovScale*100);
            end
            fprintf('\n');

            for c = 1:obj.numCams
                idx = (c-1)*6 + 1;
                pos = obj.currentChrom(idx:idx+2);
                ori = obj.currentChrom(idx+3:idx+5);
                R = eul2rotm(ori, 'XYZ');

                % Determine focal length, range, and colour per lens type
                if mod(c,2) == 0
                    f = fWide;
                    frustLen = rangeWide * fovScale;
                    col = colWide;
                    lensLabel = 'W';
                else
                    f = fNarrow;
                    frustLen = rangeNarrow * fovScale;
                    col = colNarrow;
                    lensLabel = 'N';
                end

                % Compute half-angles from pinhole model
                halfW = atan2(resolution(1)/2 * pixSize, f);
                halfH = atan2(resolution(2)/2 * pixSize, f);

                % Frustum corners in camera frame (camera looks along +Z)
                corners_cam = frustLen * [
                    tan(halfW),  tan(halfH), 1;
                   -tan(halfW),  tan(halfH), 1;
                   -tan(halfW), -tan(halfH), 1;
                    tan(halfW), -tan(halfH), 1];

                % Transform to world frame
                corners_world = (R * corners_cam')' + pos;

                % Draw frustum edges
                for e = 1:4
                    plot3([pos(1) corners_world(e,1)], ...
                          [pos(2) corners_world(e,2)], ...
                          [pos(3) corners_world(e,3)], ...
                          '-', 'Color', [col 0.5], 'LineWidth', 1);
                end

                % Draw frustum face edges
                faceOrder = [1 2 3 4 1];
                plot3(corners_world(faceOrder,1), ...
                      corners_world(faceOrder,2), ...
                      corners_world(faceOrder,3), ...
                      '-', 'Color', [col 0.6], 'LineWidth', 0.8);

                % Fill frustum face with transparency
                fill3(corners_world(1:4,1), corners_world(1:4,2), corners_world(1:4,3), ...
                    col, 'FaceAlpha', 0.08, 'EdgeColor', 'none');

                % Camera marker
                plot3(pos(1), pos(2), pos(3), 'o', 'Color', col, ...
                    'MarkerSize', 8, 'MarkerFaceColor', col);

                % Label
                text(pos(1), pos(2), pos(3)+0.3, ...
                    sprintf('C%d(%s)', c, lensLabel), ...
                    'FontSize', 9, 'FontWeight', 'bold', 'Color', col, ...
                    'HorizontalAlignment', 'center');
            end

            % Plot target points faintly
            if showTarget && isfield(obj.specs, 'Target')
                scatter3(obj.specs.Target(:,1), obj.specs.Target(:,2), obj.specs.Target(:,3), ...
                    5, [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha', 0.15);
            end

            % Draw room outline
            obj.drawRoomOutline();

            axis equal; grid on;
            xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
            title(sprintf('%d-Camera FOV (Blue=Narrow %.0fm, Red=Wide %.0fm)', ...
                obj.numCams, rangeNarrow, rangeWide));
            view(45, 30);
            legend({'', '', '', '', '', '', '', ''}, 'Visible', 'off'); % suppress auto-legend
            hold off;
        end

        %% ============================================================
        %  EXPORT SETUP SHEET
        %  ============================================================
        function exportSetupSheet(obj, filename)
            [resCost, occCost, combCost] = obj.evaluate();

            fid = fopen(filename, 'w');
            fprintf(fid, 'CAMERA SETUP SHEET\n');
            fprintf(fid, '==================\n');
            fprintf(fid, 'Generated: %s\n', char(datetime('now')));
            fprintf(fid, 'Source: %s\n', obj.sourceFile);
            fprintf(fid, 'Cameras: %d\n\n', obj.numCams);

            fprintf(fid, 'Predicted Costs:\n');
            fprintf(fid, '  Resolution Uncertainty: %.6f\n', resCost);
            fprintf(fid, '  Dynamic Occlusion:      %.6f\n', occCost);
            fprintf(fid, '  Combined:               %.6f\n\n', combCost);

            fprintf(fid, '%-5s %-6s  %-8s %-8s %-8s  %-8s %-8s %-8s  %-10s\n', ...
                'Cam', 'Lens', 'X(m)', 'Y(m)', 'Z(m)', 'Roll°', 'Pitch°', 'Yaw°', 'Wall');
            fprintf(fid, '%s\n', repmat('-', 1, 80));

            wallTol = 0.15;
            wb = obj.wallBounds;

            for c = 1:obj.numCams
                idx = (c-1)*6 + 1;
                pos = obj.currentChrom(idx:idx+2);
                ori = rad2deg(obj.currentChrom(idx+3:idx+5));

                if mod(c,2) == 0, lensStr = 'Wide'; else, lensStr = 'Narr'; end

                wallStr = 'Interior';
                if abs(pos(1) - wb(1,1)) < wallTol, wallStr = 'X-min';
                elseif abs(pos(1) - wb(1,2)) < wallTol, wallStr = 'X-max';
                elseif abs(pos(2) - wb(2,1)) < wallTol, wallStr = 'Y-min';
                elseif abs(pos(2) - wb(2,2)) < wallTol, wallStr = 'Y-max';
                end

                fprintf(fid, '%-5d %-6s  %-8.2f %-8.2f %-8.2f  %-8.1f %-8.1f %-8.1f  %-10s\n', ...
                    c, lensStr, pos(1), pos(2), pos(3), ori(1), ori(2), ori(3), wallStr);
            end

            fprintf(fid, '\nNotes:\n');
            fprintf(fid, '- Odd cameras = narrow lens (f=%.1fmm), even = wide (f=%.1fmm)\n', ...
                obj.specs.Focal*1000, obj.specs.FocalWide*1000);
            fprintf(fid, '- Positions in metres relative to room origin\n');
            fprintf(fid, '- Orientations in degrees (Roll, Pitch, Yaw — XYZ Euler)\n');
            fprintf(fid, '- Use OptiTrack Motive reprojection overlay to fine-align\n');
            fclose(fid);
            fprintf('Setup sheet exported: %s\n', filename);
        end

        %% ============================================================
        %  RESET TO ORIGINAL
        %  ============================================================
        function reset(obj)
            obj.currentChrom = obj.originalChrom;
            fprintf('Reset to original GA solution.\n\n');
            obj.printCameras();
            obj.evaluate();
        end
    end

    %% ================================================================
    %  PRIVATE METHODS
    %  ================================================================
    methods (Access = private)

        function plotCamerasInAxis(obj, chrom)
            hold on;
            for c = 1:obj.numCams
                idx = (c-1)*6 + 1;
                pos = chrom(idx:idx+2);
                ori = chrom(idx+3:idx+5);
                T = se3(eul2rotm(ori, 'XYZ'), pos);
                cam = CentralCamera('name', sprintf('cam%d', c), 'pose', T);
                cam.plot_camera('label', 'scale', 0.5);
            end
            obj.drawRoomOutline();
            axis equal; grid on;
            xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
            view(45, 30);
            hold off;
        end

        function drawRoomOutline(obj)
            % Draw capture volume (solid lines)
            vol = obj.volume;
            x = [vol(1,1) vol(1,2)];
            y = [vol(2,1) vol(2,2)];
            z = [vol(3,1) vol(3,2)];

            % Floor
            plot3([x(1) x(2) x(2) x(1) x(1)], ...
                  [y(1) y(1) y(2) y(2) y(1)], ...
                  [z(1) z(1) z(1) z(1) z(1)], ...
                  'k-', 'LineWidth', 1.0);
            % Ceiling
            plot3([x(1) x(2) x(2) x(1) x(1)], ...
                  [y(1) y(1) y(2) y(2) y(1)], ...
                  [z(2) z(2) z(2) z(2) z(2)], ...
                  'k-', 'LineWidth', 1.0);
            % Verticals
            for xi = x
                for yi = y
                    plot3([xi xi], [yi yi], z, 'k-', 'LineWidth', 0.5);
                end
            end

            % Draw room walls (dashed lines)
            wb = obj.wallBounds;
            xw = [wb(1,1) wb(1,2)];
            yw = [wb(2,1) wb(2,2)];
            zw = [wb(3,1) wb(3,2)];

            plot3([xw(1) xw(2) xw(2) xw(1) xw(1)], ...
                  [yw(1) yw(1) yw(2) yw(2) yw(1)], ...
                  [zw(1) zw(1) zw(1) zw(1) zw(1)], ...
                  'k--', 'LineWidth', 0.6, 'Color', [0.5 0.5 0.5]);
            plot3([xw(1) xw(2) xw(2) xw(1) xw(1)], ...
                  [yw(1) yw(1) yw(2) yw(2) yw(1)], ...
                  [zw(2) zw(2) zw(2) zw(2) zw(2)], ...
                  'k--', 'LineWidth', 0.6, 'Color', [0.5 0.5 0.5]);
            for xi = xw
                for yi = yw
                    plot3([xi xi], [yi yi], zw, '--', 'LineWidth', 0.4, 'Color', [0.5 0.5 0.5]);
                end
            end
        end
    end
end