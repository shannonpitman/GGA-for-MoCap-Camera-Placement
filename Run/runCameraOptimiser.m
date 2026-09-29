%% Camera Placement Optimiser for an optical MoCape tracking system
% This guided genetic algorithm optimises the arrangement of multi-camera 
% network. 
% 
% USAGE:
% The user must specify inputs below (look for any code starting with set).
% The cost function is a weighted combination of resolution uncertainty and
% dynamic occlusion handling and can also be adjusted by the user.

clc; %clear screen
clear; % clear workspace
close all;

% Add every code subfolder (GA_Core, CostFunctions, Geometry, Setup,
% Plotting, Analysis) to the MATLAB path for this session.
addProjectPaths();

%% User Inputs
% RUN CONFIGURATION: workspace, mountable regions, hardware, weights and GA
% settings all live in Setup/runConfig.m. Pick a preset and override
% anything run-specific here, e.g.
%   cfg = runConfig('optitrack_lab', 'Volume', [-3 3; -3 3; 0 3]);
%   cfg = runConfig('lowcost_tripod');
cfg = runConfig('optitrack_lab');

% CAMERA NETWORK
numCams = 7; 

% TARGET SPACE MODALITY
% 1= UAV: entire space (full flight volume)
% 2= UGV: focus on floor plane (small slab above floor; height in cfg.UGV_MaxHeight)
targetType = 2;

% Discretisation method of grid
% 1= Uniform grid (evenly spaced volume)
% 2= Normalised grid (concentrated discretisation in the centre)
targetMode = 1;

% Grid spacing [m] - x-y (in-plane). For UGV the z spacing is cfg.UGV_ZSpacing.
spacing = 1;

% COST FUNCTION:
% 1 = Resolution Uncertainty only
% 2 = Dynamic Occlusion only  
% 3 = Combined (weighted by cfg.Weights)
costFunctionType = 3;

% WARM-START
% To warm-start the GA from a previously found chromosome, set
% warmStartUsed = true and assign warmStartBestSol to a saved
% chromosome (1 x 6*numCams row vector — same convention as
% saveData.BestSolution.Chromosome). The chromosome MUST come from
% a run with the same cost functions; older chromosomes were optimised
% against a different cost surface and will mislead the GA.
warmStartUsed    = false;
warmStartBestSol = [];      % e.g. load(...).saveData.BestSolution.Chromosome

%% Set-up
[specs, problem, params] = buildRunSpecs(cfg, numCams, costFunctionType, ...
    targetType, targetMode, spacing);

if warmStartUsed
    if isempty(warmStartBestSol) || numel(warmStartBestSol) ~= problem.nVar
        error('runCameraOptimiser:BadWarmStart', ...
            'warmStartUsed=true requires warmStartBestSol of length %d.', problem.nVar);
    end
    perturbed = Mutate(warmStartBestSol, 1, 0.5);                 % perturb all genes
    perturbed = min(max(perturbed, problem.VarMin), problem.VarMax);
    specs.warmStart = true;
    specs.warmChromosomes = [warmStartBestSol; perturbed];
end

%% Run GA
tic; % start timer
out = RunGA(problem, params, specs);
elapsedTime = toc; % end timer

fprintf('\n :> Optimisation Complete :>\n');
fprintf('Best Cost: %.6f\n', out.bestsol.Cost);
fprintf('Computation Time: %.2f min\n', elapsedTime / 60);

%% Results 
coverageStats = visualizeCameraCoverage(out, specs);
plotResults(out, specs, params, elapsedTime);
saveResults(out, specs, params, elapsedTime, costFunctionType, warmStartUsed, coverageStats)