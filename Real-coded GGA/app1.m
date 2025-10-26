clc; %clear screen 
clear; % clear workspace
close all;

%% Design Specifications
specs.Cams = 6; %Number of Cameras
specs.Resolution = [640 480]; %VGA resolution
specs.PixelSize = 1.4e-6; %Square Pixel Size
specs.PrincipalPoint = [specs.Resolution(1)/2, specs.Resolution(2)/2];
specs.Focal = 0.0028; %focal length [m]

% Target space is an uniformly discretised grid within the flight volume 
% This workspace volume matches the available dimensions of the MS.G flight envelope 
flight_envelope = [-4 4; -4 4; 0 4.5]; %m
spacing = 0.5;
x_marker = flight_envelope(1,1):spacing:flight_envelope(1,2);
y_marker = flight_envelope(2,1):spacing:flight_envelope(2,2);
z_marker = flight_envelope(3,1):spacing:flight_envelope(3,2);
[X,Y,Z] = meshgrid(x_marker, y_marker, z_marker);

specs.Target = [X(:), Y(:), Z(:)];

numSections = floor(specs.Cams/2);
nx = ceil(sqrt(numSections)); %number of divisions in x direction of grid 
ny = ceil(numSections /nx); %number of divisions in y direction of grid 
x_div = linspace(flight_envelope(1,1), flight_envelope(1,2), nx+1);
y_div = linspace(flight_envelope(2,1), flight_envelope(2,2), ny+1);
z_mid = mean(flight_envelope(3,:)); %mid-height
count = 1;
section_centres = zeros(numSections,3);
for i = 1:length(x_div)-1
    for j = 1:length(y_div)-1
        if count <= numSections
            cx = mean([x_div(i), x_div(i+1)]);
            cy = mean([y_div(j), y_div(j+1)]);
            cz = z_mid;
            section_centres(count,:) = [cx, cy, cz];
            count = count+1;
        end
    end
end

specs.SectionCentres = section_centres;

%% Problem Definition

problem.CostFunction = @resUncertainty; % Objective function 
%Bounds for the 6 design variables (search space) [Xc, Yc, Zc, alpha, beta, gamma]
cameraLowerBounds = [-5 -4.5  0   -pi -pi -pi];
cameraUpperBounds = [ 5  4.5  4.8  pi  pi  pi];
problem.VarMin = repmat(cameraLowerBounds,1,specs.Cams);
problem.VarMax = repmat(cameraUpperBounds,1,specs.Cams);
problem.nVar = 6* specs.Cams;

%% GA Parameters

params.MaxIt = 60;
params.nPop = 300;

params.beta = 1;
params.pC = 1;
params.gamma = 0.1;
params.mu = 0.02; %probability of mutation
params.sigma = 00.1;

%% Run GA
out = RunGA(problem, params, specs);

%% Results 

figure;
%plot(out.bestcost, 'LineWidth',2);
semilogy(out.bestcost, 'LineWidth',2); % y axis has logarithmic scale and x-axis is linear 
xlabel('Iterations');
ylabel('Best Cost');
grid on;