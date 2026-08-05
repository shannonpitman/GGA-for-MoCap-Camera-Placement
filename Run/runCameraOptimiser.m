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
% CAMERA NETWORK
numCams = 7; 

% WORKSPACE VOLUME [m] & MOUNTING CONSTRAINTS
volume = [-4 4; -4 4; 0 4]; %Matches the available dimensions of the MS.G flight envelope 
%Bounds for the 6 design variables (search space) [Xc, Yc, Zc, alpha, beta, gamma]
cameraLowerBounds = [-5 -4.5 0  -pi -pi/2 -pi]; %wall-mounting constraints lower bounds
cameraUpperBounds = [ 5  4.5 4.8 pi  pi/2  pi]; %wall-mounting constraints upper bounds

% TARGET SPACE MODALITY
% 1= UAV: entire space (full flight volume)
% 2= UGV: focus on floor plane (small slab above floor)
targetType = 2;

% Discretisation method of grid
% 1= Uniform grid (evenly spaced volume)
% 2= Normalised grid (concentrated discretisation in the centre)
targetMode = 1;

% Grid spacing [m] - x-y (in-plane). For UGV the z spacing is decoupled
% so that x-y can be coarsened without losing slab layers.
spacing = 1;

% UGV max expected height of markers on UGV [m]
UGV_maxHeight = 0.5;     % only used if targetType == 2
UGV_zSpacing  = 0.25;    % z step on the UGV slab; 0.25 m -> 3 layers

if targetType == 2
    volume(3,:) = [0, UGV_maxHeight];
    % Anisotropic spacing: keep x-y as supplied, fix z to UGV_zSpacing.
    targetSpacing = [spacing, spacing, min(UGV_zSpacing, UGV_maxHeight)];
else
    targetSpacing = spacing;   % isotropic for UAV
end

% COST FUNCTION:
% 1 = Resolution Uncertainty only
% 2 = Dynamic Occlusion only  
% 3 = Combined (weighted)

costFunctionType = 3;
problem = setupProblem(numCams, costFunctionType, cameraUpperBounds, cameraLowerBounds);

% GA PARAMETERS
% Population scales with chromosome length: nPop = (params per camera) *
% numCams * 10 = problem.nVar * 10. numParams is read from the per-camera
% bound vector so it stays correct if the camera model changes.
numParams      = numel(cameraLowerBounds);        % design variables per camera (6)
maxGenerations = 100;
populationSize = numCams * numParams * 10;         % e.g. 7 cams -> 420

% WARM-START
% To warm-start the GA from a previously found chromosome, set
% warmStartUsed = true and assign warmStartBestSol to a saved
% chromosome (1 x 6*numCams row vector — same convention as
% saveData.BestSolution.Chromosome). The chromosome MUST come from
% a post-bugfix run; pre-fix chromosomes were optimised against the
% corrupted cost surface and will mislead the GA.
warmStartUsed    = false;
warmStartBestSol = [];      % e.g. load(...).saveData.BestSolution.Chromosome

if warmStartUsed
    if isempty(warmStartBestSol) || numel(warmStartBestSol) ~= 6*numCams
        error('runCameraOptimiser:BadWarmStart', ...
            'warmStartUsed=true requires warmStartBestSol of length %d.', ...
            6*numCams);
    end
    perturbed = Mutate(warmStartBestSol, 1, 0.5);                 % perturb all genes
    perturbed = max(perturbed, problem.VarMin);
    perturbed = min(perturbed, problem.VarMax);
    warmStartChromosome = [warmStartBestSol; perturbed];
else
    warmStartChromosome = [];
end

%% Set-up
% Hardware (input camera intrinsics)
specs = setupHardwareSpecs(numCams);

specs.warmStart = warmStartUsed;
specs.warmChromosomes = warmStartChromosome;

%Weights for combined cost function (only used if costFunctionType = 3)
%Weights need to sum to 1
specs.WeightUncertainty = 0.5; % Resolution uncertainty weight
specs.WeightOcclusion = 0.5; % Dynamic occlusion weight

% Robotics Modality
specs.TargetType = targetType;

% Target Space
specs.TargetMode = targetMode;
specs.Target = generateTargetSpace(volume, targetMode, targetSpacing);
specs.NumPoints = size(specs.Target,1);
specs.spacing = spacing;          % scalar x-y; kept for log
if targetType == 2
    specs.spacingZ = UGV_zSpacing;
end

% Section centres for guided initial population
specs.SectionCentres = generateSectionCentres(numCams, volume);

% Cost function parameters
specs = setupCostParams(specs);

% GA Parameters
params = setupGAparams(maxGenerations, populationSize);

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