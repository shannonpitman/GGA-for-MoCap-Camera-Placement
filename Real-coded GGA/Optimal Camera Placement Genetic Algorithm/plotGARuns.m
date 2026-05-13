%% Generates figures from GA results from batch-run
clear; clc; close all;

% Make sure every code subfolder is on the MATLAB path.
projectRoot = fileparts(mfilename('fullpath'));
addProjectPaths();

%% File Management (post-restructure layout)
logFile   = fullfile(projectRoot, 'Results', 'Logs', 'GGA_RunsLog.mat'); % Master run log
runDir    = fullfile(projectRoot, 'Results');   % Parent of <N>Cams subfolders
outputDir = fullfile(projectRoot, 'figures');   % Where to save PDFs

if ~isfolder(outputDir)
    mkdir(outputDir);
end


%% Experimental conditions
% Camera counts tested
cameraRange = [6, 7, 8];
 
% Cost functions: 1=ResUnc, 2=DynOcc, 3=Combined
costFuncs   = [1, 2, 3];
cfLabels    = {'CF1_ResUnc', 'CF2_DynOcc', 'CF3_Combined'};
 
% Target types: 1=UAV, 2=UGV
targetTypes = [1, 2];
ttLabels    = {'UAV', 'UGV'};
 
% Grid modes: 1=Uniform, 2=Normal
gridModes   = [1, 2];
gmLabels    = {'GM1_Uniform', 'GM2_Normal'};
 
% Spacings tested (metres)
spacings    = [0.50; 1; 1.50; 2];
spLabels    = arrayfun(@(s) sprintf('sp%.0fcm', s*100), spacings, 'Uni', false);
 
% Representative camera count for convergence curves
repCams = 7;
 
% Common filter arguments used across most plots
commonArgs = {'LogFile', logFile};
 
 
%% Figure Set 1: Convergence Curves
%  One per (CostFunction x TargetType x GridMode x Spacing).
%  Each figure shows repeated runs under identical conditions.
%  This is the key methodology evidence that the GA converges reliably.
 
fprintf('\n===== CONVERGENCE CURVES =====\n');
fprintf('  (one per condition, %dC representative)\n\n', repCams);
 
for cf = costFuncs
    for tt = targetTypes
        for gm = gridModes
            for si = 1:length(spacings)
                sp = spacings(si);
                tag = sprintf('convergence_%s_%s_%s_%s_%dC', ...
                    cfLabels{cf}, ttLabels{tt}, gmLabels{gm}, spLabels{si}, repCams);
                try
                    plotGA_Convergence(commonArgs{:}, ...
                        'CostFunction', cf, ...
                        'TargetType',   tt, ...
                        'GridMode',     gm, ...
                        'Spacing',      sp, ...
                        'NumCameras',   repCams, ...
                        'RunDir',       runDir, ...
                        'ShowTopTen',   true, ...
                        'ShowAvgCost',  false, ...
                        'LogScale',     true, ...
                        'SaveAs', fullfile(outputDir, tag));
                catch ME
                    fprintf('  Skipped %s: %s\n', tag, ME.message);
                end
            end
        end
    end
end
 
 
%% ========== Figure Set 2: Cost Box Plots (by Camera Count) ==========
%  One figure per (CostFunction x GridMode x Spacing), split by TargetType.
%  Shows how cost decreases with more cameras, UAV vs UGV side by side.
 
fprintf('\n===== COST BOX PLOTS (split by Target Type) =====\n');
 
for cf = costFuncs
    for gm = gridModes
        for si = 1:length(spacings)
            sp = spacings(si);
            tag = sprintf('cost_byTT_%s_%s_%s', ...
                cfLabels{cf}, gmLabels{gm}, spLabels{si});
            try
                plotGA_CostBoxPlots(commonArgs{:}, ...
                    'CostFunction', cf, ...
                    'GridMode',     gm, ...
                    'Spacing',      sp, ...
                    'SplitBy',      'TargetType', ...
                    'SaveAs', fullfile(outputDir, tag));
            catch ME
                fprintf('  Skipped %s: %s\n', tag, ME.message);
            end
        end
    end
end
 
 
%% ========== Figure Set 3: Warm-Start vs Cold-Start ==========
%  One figure PER CAMERA COUNT is emitted by plotGA_WarmColdEffect
%  internally, so the outer loop only enumerates the experimental
%  conditions that distinguish the warm-start *story*. Restricted to
%  GM=1 (Uniform): GM=2 (Normal) is retained only for coverage /
%  best-configuration interest and would clutter the warm-cold result
%  set without adding evidence for the warm-start claim.

fprintf('\n===== WARM-START EFFECT =====\n');

warmColdGM = 1;  % only Uniform grid for this figure set

for cf = costFuncs
    for tt = targetTypes
        for si = 1:length(spacings)
            sp = spacings(si);
            tag = sprintf('warmcold_%s_%s_%s_%s', ...
                cfLabels{cf}, ttLabels{tt}, gmLabels{warmColdGM}, spLabels{si});
            try
                plotGA_WarmColdEffect(commonArgs{:}, ...
                    'CostFunction', cf, ...
                    'TargetType',   tt, ...
                    'GridMode',     warmColdGM, ...
                    'Spacing',      sp, ...
                    'SaveAs', fullfile(outputDir, tag));
            catch ME
                fprintf('  Skipped %s: %s\n', tag, ME.message);
            end
        end
    end
end
 
 
%% ========== Figure Set 4: Computation Time ==========
%  Split by target type (UAV full volume vs UGV floor slab).
%  Controlled per (GridMode x Spacing) since point count drives time.
 
fprintf('\n===== COMPUTATION TIME =====\n');
 
for gm = gridModes
    for si = 1:length(spacings)
        sp = spacings(si);
 
        % Split by target type
        tag = sprintf('comptime_byTT_%s_%s', gmLabels{gm}, spLabels{si});
        try
            plotGA_ComputationTime(commonArgs{:}, ...
                'GridMode',  gm, ...
                'Spacing',   sp, ...
                'SplitBy',   'TargetType', ...
                'SaveAs', fullfile(outputDir, tag));
        catch ME
            fprintf('  Skipped %s: %s\n', tag, ME.message);
        end
 
        % Split by cost function (shows relative cost of each evaluation)
        tag = sprintf('comptime_byCF_%s_%s', gmLabels{gm}, spLabels{si});
        try
            plotGA_ComputationTime(commonArgs{:}, ...
                'GridMode',  gm, ...
                'Spacing',   sp, ...
                'SplitBy',   'CostFunction', ...
                'SaveAs', fullfile(outputDir, tag));
        catch ME
            fprintf('  Skipped %s: %s\n', tag, ME.message);
        end
    end
end
 
 
%% ========== Figure Set 5: Factor Effects (Grid Mode comparison) ==========
%  Side-by-side: TargetType effect (left) + GridMode effect (right).
%  One per (CostFunction x Spacing). This is the ONE figure set that
%  deliberately compares across grid modes, so spacing must still match.
 
fprintf('\n===== FACTOR EFFECTS =====\n');
 
for cf = costFuncs
    for si = 1:length(spacings)
        sp = spacings(si);
        tag = sprintf('factors_%s_%s', cfLabels{cf}, spLabels{si});
        try
            plotGA_FactorEffects(commonArgs{:}, ...
                'CostFunction', cf, ...
                'Spacing',      sp, ...
                'SaveAs', fullfile(outputDir, tag));
        catch ME
            fprintf('  Skipped %s: %s\n', tag, ME.message);
        end
    end
end
 
 
%% ========== Summary ==========
fprintf('\n===== DONE =====\n');
pdfFiles = dir(fullfile(outputDir, '*.pdf'));
fprintf('Generated %d PDF figures in %s/\n', length(pdfFiles), outputDir);
for i = 1:length(pdfFiles)
    fprintf('  %s\n', pdfFiles(i).name);
end
 
fprintf('\nFigure naming convention:\n');
fprintf('  <type>_<CF>_<TT>_<GM>_<spacing>[_<cameras>].pdf\n');
fprintf('  e.g. convergence_CF3_Combined_UAV_GM1_Uniform_sp50cm_7C.pdf\n');
fprintf('\nInclude in LaTeX via:\n');
fprintf('  \\includegraphics[width=\\textwidth]{figures/<filename>.pdf}\n');
 