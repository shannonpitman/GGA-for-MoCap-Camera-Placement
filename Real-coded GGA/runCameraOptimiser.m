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
targetType = 2;

% Discretisation method of grid
% 1= Uniform grid (evenly spaced volume)
% 2= Normalised grid (concentrated discretisation in the centre)
targetMode = 1;

% Grid spacing [m]
spacing = 0.5; %increase for faster evaluation, increase for UGV mode

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
problem = setupProblem(numCams, costFunctionType, cameraUpperBounds, cameraLowerBounds);

% GA PARAMETERS
maxGenerations = 150;
populationSize = 100;

% WARM-START
warmStartUsed = false; % Set this to true when you use warm-start
warmStartBestSol = [-3.94453694866477	2.80270658457369	4.70420099523774	3.14028385018124	0.865374203482423	1.52808794739770	4.81414500579010	1.49344760628434	4.07726936557076	3.13638271864850	-0.421481339363494	-2.45179618042072	-4.95492972522762	-2.48343183175008	4.54438986956409	-1.91517534692494	0.350745443390892	-1.63869078772669	-4.87109356019903	1.97091349231205	3.73307831294323	3.07387024977452	0.275222610101504	2.32742052217727	1.98977952260285	-4.50000000000000	4.24773888189726	-2.08157684817404	-0.170059426054804	0.423818681832194	4.74229673914245	1.50991375796534	4.62881140632982	-3.14159265358979	-0.227255198303810	3.14159265358979	-0.484212964041819	-4.47751014807785	2.13048050788250	-1.72318139201811	-0.0469945210243474	-1.1266562383236155]; % Uncomment and insert chromosome to seed pop(1)
perturbed = Mutate(warmStartBestSol, 1, 0.5); %mutates all genes 
perturbed = max(perturbed, problem.VarMin);
perturbed = min(perturbed, problem.VarMax);
warmStartChromosome = [warmStartBestSol; perturbed];
>>>>>>> 3c2db45fdf89886dd088769d5b0ce831e4cb2f08

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
specs.Target = generateTargetSpace(volume, targetMode, spacing);
specs.NumPoints = size(specs.Target,1);

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