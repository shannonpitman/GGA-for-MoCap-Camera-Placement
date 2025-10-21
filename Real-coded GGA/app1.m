clc; %clear screen 
clear; % clear workspace
close all;

%% Design Specifications
specs.Cams = 4; %Number of Cameras
% Target space is an uniformly discretised grid within the flight volume 
% This workspace volume matches the available dimensions of the MS.G flight envelope 
flight_envelope = [-4 4; -4 4; 0 4.5]; %m
x_marker = flight_envelope(1,1):spacing:flight_envelope(1,2);
y_marker = flight_envelope(2,1):spacing:flight_envelope(2,2);
z_marker = flight_envelope(3,1):spacing:flight_envelope(3,2);
[X,Y,Z] = meshgrid(x_marker, y_marker, z_marker);
specs.Target = [X(:), Y(:), Z(:)];

%% Problem Definition

problem.CostFunction = @resUncertainty; % Objective function 
problem.nVar = 6; % Amount of design variables (search space) [Xc, Yc, Zc, alpha, beta, gamma]
problem.VarMin = [-5 -4.5  0   -pi -pi -pi]; %lower bounds of variables 
problem.VarMax = [ 5  4.5  4.8  pi  pi  pi]; % upper bounds of variables

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