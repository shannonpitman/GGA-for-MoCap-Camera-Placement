clc; %clear screen 
clear; % clear workspace
close all;

%% Design Specifications
specs.Cams = 7; %Number of Cameras
specs.Resolution = [640 480]; %VGA resolution
specs.PixelSize = 1.4e-6; %Square Pixel Size
specs.PrincipalPoint = [specs.Resolution(1)/2, specs.Resolution(2)/2];
specs.Focal = 0.0028; %focal length [m]

% Target space is an uniformly discretised grid within the flight volume 
% This workspace volume matches the available dimensions of the MS.G flight envelope 
flight_envelope = [-3 3; -3 3; 0 2]; %m -> desired space: [-4 4; -4 4; 0 4.5]
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
            section_centres(count,:) = [(x_div(i) + x_div(i+1))/2, (y_div(j)+ y_div(j+1))/2, z_mid];
            count = count+1;
        end
    end
end

specs.SectionCentres = section_centres;

%% Problem Definition
% Choose cost function:
% 1 = Resolution Uncertainty only
% 2 = Dynamic Occlusion only  
% 3 = Combined (weighted)

costFunctionType = 1; % Change this to select cost function

%Weights for combined cost function (only used if costFunctionType = 3)
%Weights need to sum to 1
specs.WeightUncertainty = 0.7; % Resolution uncertainty weight
specs.WeightOcclusion = 0.3; % Dynamic occlusion weight

% Select cost function
switch costFunctionType
    case 1
        problem.CostFunction = @resUncertainty;
    case 2
        problem.CostFunction = @dynamicOcclusion;
    case 3
        problem.CostFunction = @combinedCostFunction;
end

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
params.mu = 0.05; %probability of mutation
params.sigma = 00.1;

%% Run GA
tic; % start timer
out = RunGA(problem, params, specs);
elapsedTime = toc; % end timer

%% Results 
currentDateTime = datetime('now');
dateTimeStr = string(currentDateTime, 'yyyyMMdd_HHmmSS');
filename = sprintf('%dCams_Run_%s.mat', specs.Cams, dateTimeStr);

% Prepare data to save
saveData.BestSolution = out.bestsol;
saveData.BestCost = out.bestsol.Cost;
saveData.CameraConfiguration = reshape(out.bestsol.Chromosome, 6, specs.Cams)';
saveData.Specifications = specs;
saveData.GAParams = params;
saveData.ConvergenceHistory = out.bestcost;
saveData.ElapsedTime = elapsedTime;
saveData.Timestamp = currentDateTime;

% Save camera positions and orientations (best solution configuaratrion)
for i = 1:specs.Cams
    chromStartIdx = (i-1)*6+1;
    chromEndIdx = i*6;
    saveData.Cameras(i).Position = out.bestsol.Chromosome(chromStartIdx:chromStartIdx+2);
    saveData.Cameras(i).Orientation = out.bestsol.Chromosome(chromEndIdx-2:chromEndIdx);
    saveData.Cameras(i).OrientationDegrees = rad2deg(saveData.Cameras(i).Orientation);
end
% Save to MAT file
save(filename, 'saveData');
fprintf('\n Optimization Complete :>\n');
fprintf('Results saved to: %s\n', filename);
fprintf('Best Cost: %.6f\n', saveData.BestCost);
fprintf('Computation Time: %.2f seconds\n', elapsedTime);

% Also save a human-readable text summary
txtFilename = sprintf('%dCams_Run_%s.txt', specs.Cams, dateTimeStr);
fid = fopen(txtFilename, 'w');
fprintf(fid, 'Genetic Algorithm Camera Placement Results\n');
fprintf(fid, '==========================================\n\n');
fprintf(fid, 'Timestamp: %s\n', char(currentDateTime));
fprintf(fid, 'Number of Cameras: %d\n', specs.Cams);
fprintf(fid, 'Best Cost (Uncertainty): %.6f\n', saveData.BestCost);
fprintf(fid, 'Computation Time: %.2f seconds\n\n', elapsedTime);

fprintf(fid, 'Camera Configurations:\n');
fprintf(fid, '---------------------\n');
for i = 1:specs.Cams
    fprintf(fid, '\nCamera %d:\n', i);
    fprintf(fid, '  Position (m): [%.3f, %.3f, %.3f]\n', saveData.Cameras(i).Position);
    fprintf(fid, '  Orientation (rad): [%.3f, %.3f, %.3f]\n', saveData.Cameras(i).Orientation);
    fprintf(fid, '  Orientation (deg): [%.1f, %.1f, %.1f]\n', saveData.Cameras(i).OrientationDegrees);
end

fclose(fid);
fprintf('Summary saved to: %s\n', txtFilename);


figure('Name', 'GA Convergence', 'Position', [100, 100, 800, 600]);
semilogy(out.bestcost, 'LineWidth', 2); % y axis has logarithmic scale and x-axis is linear 
xlabel('Iterations');
ylabel('Best Cost (logarithmic)');
title(sprintf('GA Convergence - %d Cameras', specs.Cams))
grid on;
% convergence statistics
text(0.6*params.MaxIt, max(out.bestcost)*0.5, sprintf('Final Cost: %.4f\nTime: %.1fs', saveData.BestCost, elapsedTime),'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k');

plotFilename = sprintf('%dCams_Run_%s_convergence.png', specs.Cams, dateTimeStr);
saveas(gcf, plotFilename);
fprintf('Convergence plot saved to: %s\n', plotFilename);

% Camera Visualisation
cameras = cell(specs.Cams, 1);
for i = 1:specs.Cams
    T = se3(eul2rotm(saveData.Cameras(i).Orientation, "XYZ"), saveData.Cameras(i).Position);
    cameras{i} = CentralCamera(name="cam"+i, pose=T);
end

figure('Name', 'Camera Configuration', 'Position', [950, 100, 800, 600]);
hold on

% Plot cameras
for i = 1:specs.Cams
    cameras{i}.plot_camera('label', scale = 0.5)
end

axis('equal')
grid('on')
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
title(sprintf('Optimal Camera Placement - Cost: %.4f', saveData.BestCost))
view(45, 30)
hold off

cameraPlotFilename = sprintf('%dCams_Run_%s_cameras.png', specs.Cams, dateTimeStr);
saveas(gcf, cameraPlotFilename); %saves as png

%Animation Plot 
% Create animation of evolving camera configurations
animateFig = figure('Name', 'Camera Configuration Evolution');
for frameIdx = params.MaxIt
    hold on;
    chromosome = out.bestChromosomes(frameIdx, :);
    % Update camera positions for the next frame in the animation
    for i = 1:specs.Cams
        chromStartIdx = (i-1)*6+1;
        chromEndIdx = i*6;
        T = se3(eul2rotm(chromosome(chromEndIdx-2:chromEndIdx), "XYZ"), chromosome(chromStartIdx:chromStartIdx+2));
        cameras{i} = CentralCamera(name="cam"+i, pose=T);
        cameras{i}.plot_camera('label', scale = 0.5);
    end
    drawnow;
    frame = getframe(animateFig);
    im{frameIdx} = frame2im(frame);
end
%save animation as a gif
% filenameAnimate = sprintf('%dEvolutionofCams_Run_%s.gif', specs.Cams, dateTimeStr);
% for idx  = 1:params.MaxIt 
%     [A,map] = rgb2ind(im{idx}, 256);
%     if idx == 1
%         imwrite(A,map,filenameAnimate, "gif", LoopCount= inf, DelayTime=1);
%     else
%         imwrite(A, map, filenameAnimate, "gif", 'WriteMode', 'append', 'DelayTime', 1);
%     end
% end
% fprintf('Animation saved to: %s\n', filenameAnimate);

%% Target Space Visualisation (Coverage)
% out.bestsol.Chromosome = [-0.150, 3.951, 4.800, 2.489, -0.587,-1.333,...
%                             5.000, 4.285, 1.453, 1.336, -0.873, 0.293, ...
%                             -2.644, 4.500, 3.196, 1.876, 0.486, 0.393, ...
%                             -2.578, -4.500, 3.623, -1.961, -0.048, 0.029, ...
%                             5.000, 4.418, 1.525, 1.299, -0.721, -0.398, ...
%                             -4.516, 4.500, 1.599, 1.436, 0.763, 0.645, ...
%                             -2.751, -4.500, 2.214, -1.556, 0.484, -0.117];
% currentDateTime = datetime('now');
% dateTimeStr = string(currentDateTime, 'yyyyMMdd_HHmmSS');
% saveData.BestCost = 73.73;
figTitle = sprintf('Camera Coverage - %d Cameras (Cost: %.4f)', specs.Cams, saveData.BestCost);
visualizeCameraCoverage(out.bestsol.Chromosome, specs, figTitle);

coveragePlotFilename = sprintf('%dCams_Run_%s_coverage.png', specs.Cams, dateTimeStr);
saveas(gcf, coveragePlotFilename);
fprintf('Coverage plot saved to: %s\n', coveragePlotFilename);
