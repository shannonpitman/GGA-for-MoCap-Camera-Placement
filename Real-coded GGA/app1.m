clc; %clear screen 
clear; % clear workspace
close all;
%% Specify Path

onedrive = fullfile(getenv('USERPROFILE'), 'OneDrive - University of Cape Town');
addpath(fullfile(onedrive, 'MSc', 'MATLAB', 'RVC3-MATLAB'));
addpath(fullfile(onedrive, 'MSc', 'MATLAB', 'Simulation', 'Optimising Camera Placement'));
%% Number of Cameras
numCams = 6; 
%% Problem Definition

problem.CostFunction = @resUncertainty; % Objective function 
problem.nVar = 5; % Amount of design variables (search space) 
problem.VarMin = [-10 -10 -5 -1 -5]; %lower bounds of variables 
problem.VarMax = [ 10  10  5  1 8];

%% GA Parameters

params.MaxIt = 60;
params.nPop = 300;

params.beta = 1;
params.pC = 1;
params.gamma = 0.1;
params.mu = 0.02; %probability of mutation
params.sigma = 00.1;

%% Run GA
out = RunGA(problem, params);

%% Results 

figure;
%plot(out.bestcost, 'LineWidth',2);
semilogy(out.bestcost, 'LineWidth',2); % y axis has logarithmic scale and x-axis is linear 
xlabel('Iterations');
ylabel('Best Cost');
grid on;