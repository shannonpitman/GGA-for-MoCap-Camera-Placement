%% analyseBest7CamConfigs  Driver for the 7-camera best-config analysis suite.
%
% For each scenario (UAV, UGV) under CF3 / Uniform grid / 7 cameras, this
% script runs the four analysis figures introduced 2026-05-13:
%
%   1. plotHeatmap_GAvsOptiTrack       — coverage (count) heat-map
%   2. plotCostField_GAvsOptiTrack     — per-point uncertainty + occlusion
%   3. plotConfigFOV_GAvsOptiTrack     — camera-pose + FOV frustum
%   4. plotBaselineAngles_GAvsOptiTrack — pairwise triangulation-angle hist
%   5. plotAngleSensitivity             — roll/pitch/yaw sweep on one camera
%
% Outputs land in <project>/figures/ as vector PDFs.

clear; clc; close all;

projectRoot = fileparts(mfilename('fullpath'));
addProjectPaths();

logFile   = fullfile(projectRoot, 'Results', 'Logs', 'GGA_RunsLog.mat');
outputDir = fullfile(projectRoot, 'figures');
if ~isfolder(outputDir)
    mkdir(outputDir);
end

% Conditions shared by every figure in this batch
sharedArgs = { ...
    'LogFile',      logFile, ...
    'GridMode',     1, ...          % Uniform
    'Spacing',      1.0, ...        % 1 m
    'CostFunction', 3, ...          % Combined
    'NumCameras',   7};

targetTypes = [1, 2];           % UAV, UGV
ttNames     = {'UAV', 'UGV'};

% Choose which camera index to sweep for the angle-sensitivity plot.
% Default 1, but you can override per scenario if a particular camera
% looked more interesting in the FOV plot.
sweepCamPerScenario = [1, 1];

for k = 1:length(targetTypes)
    tt   = targetTypes(k);
    name = ttNames{k};
    fprintf('\n===== %s =====\n', name);

    % 1) Coverage count heat-map
    try
        plotHeatmap_GAvsOptiTrack(sharedArgs{:}, 'TargetType', tt, ...
            'SaveAs', fullfile(outputDir, sprintf('CoverageHeatmap_GAvsOpti_%s', name)));
    catch ME
        fprintf('  plotHeatmap_GAvsOptiTrack failed: %s\n', ME.message);
    end

    % 2) Per-point uncertainty + occlusion field
    try
        plotCostField_GAvsOptiTrack(sharedArgs{:}, 'TargetType', tt, ...
            'SaveAs', fullfile(outputDir, sprintf('CostField_GAvsOpti_%s', name)));
    catch ME
        fprintf('  plotCostField_GAvsOptiTrack failed: %s\n', ME.message);
    end

    % 3) Camera pose + FOV frustum
    try
        plotConfigFOV_GAvsOptiTrack(sharedArgs{:}, 'TargetType', tt, ...
            'SaveAs', fullfile(outputDir, sprintf('ConfigFOV_GAvsOpti_%s', name)));
    catch ME
        fprintf('  plotConfigFOV_GAvsOptiTrack failed: %s\n', ME.message);
    end

    % 4) Pairwise baseline-angle histogram
    try
        plotBaselineAngles_GAvsOptiTrack(sharedArgs{:}, 'TargetType', tt, ...
            'SaveAs', fullfile(outputDir, sprintf('BaselineAngles_GAvsOpti_%s', name)));
    catch ME
        fprintf('  plotBaselineAngles_GAvsOptiTrack failed: %s\n', ME.message);
    end

    % 5) Single-camera angle sensitivity
    try
        plotAngleSensitivity(sharedArgs{:}, 'TargetType', tt, ...
            'CameraIndex', sweepCamPerScenario(k), ...
            'SweepRange',  45, 'SweepStep', 3, ...
            'ShowAllCams', false, ...
            'SaveAs', fullfile(outputDir, ...
                sprintf('AngleSensitivity_%s_cam%d', name, sweepCamPerScenario(k))));
    catch ME
        fprintf('  plotAngleSensitivity failed: %s\n', ME.message);
    end
end

fprintf('\n===== DONE =====\n');
pdfFiles = dir(fullfile(outputDir, '*.pdf'));
fprintf('PDFs in %s/:\n', outputDir);
for i = 1:length(pdfFiles)
    fprintf('  %s\n', pdfFiles(i).name);
end
