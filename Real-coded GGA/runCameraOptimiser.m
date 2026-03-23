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

%% User Inputs
% CAMERA NETWORK
numCams = 7; 

% WORKSPACE VOLUME [m] & MOUNTING CONSTRAINTS
volume = [-4 4; -4 4; 0 4]; %Matches the available dimensions of the MS.G flight envelope 
%Bounds for the 6 design variables (search space) [Xc, Yc, Zc, alpha, beta, gamma]
cameraLowerBounds = [-5 -4.5 0  -pi -pi -pi]; %wall-mounting constraints lower bounds
cameraUpperBounds = [ 5  4.5 4.8 pi  pi  pi]; %wall-mounting constraints upper bounds

% TARGET SPACE MODALITY
% 1= UAV: entire space (full flight volume)
% 2= UGV: focus on floor plane (small slab above floor)
targetType = 1;

% Discretisation method of grid
% 1= Uniform grid (evenly spaced volume)
% 2= Normalised grid (concentrated discretisation in the centre)
targetMode = 1;

% Grid spacing [m]
spacing = 1; %increase for faster evaluation

% UGV max expected height of markers on UGV [m]
UGV_maxHeight = 0.5; %only used if targetType =2
if targetType ==2 %UGV max height in workspace volume
    volume(3,:) = [0, UGV_maxHeight];
    if spacing > UGV_maxHeight
        spacing = UGV_maxHeight/2; %at least two layers
    end
end

% COST FUNCTION:
% 1 = Resolution Uncertainty only
% 2 = Dynamic Occlusion only  
% 3 = Combined (weighted)

costFunctionType = 3;

% GA PARAMETERS
maxGenerations = 150;
populationSize = 150;

% WARM-START
warmStartUsed = false; % Set this to true when you use warm-start
% warmStartChromose = []; % Uncomment and insert chromosome to seed pop(1)


%% Set-up
% Hardware (input camera intrinsics)
specs = setupHardwareSpecs(numCams);

%Weights for combined cost function (only used if costFunctionType = 3)
%Weights need to sum to 1
specs.WeightUncertainty = 0.5; % Resolution uncertainty weight
specs.WeightOcclusion = 0.5; % Dynamic occlusion weight

% Target Space
specs.Target = generateTargetSpace(volume, targetMode, spacing);
specs.NumPoints = size(specs.Target,1);

% Section centres for guided initial population
specs.SectionCentres = generateSectionCentres(numCams, volume);

% Cost function parameters
specs = setupCostParams(specs);
problem = setupProblem(numCams, costFunctionType, cameraUpperBounds, cameraLowerBounds);

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