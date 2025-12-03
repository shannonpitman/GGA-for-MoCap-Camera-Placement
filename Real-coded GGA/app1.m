clc; %clear screen 
clear; % clear workspace
close all;

%% Specifications
specs.Cams = 7; %Number of Cameras
specs.Resolution = [640 480]; %VGA resolution
specs.npix = [640 480]; %VGA resolution
specs.PixelSize = 1.4e-6; %Square Pixel Size
specs.PrincipalPoint = [specs.Resolution(1)/2, specs.Resolution(2)/2];
specs.Focal = 0.0028; %focal length [m]

% Target space is a uniformly discretised grid within the flight volume 
% This workspace volume matches the available dimensions of the MS.G flight envelope 
flight_envelope = [-3 3; -3 3; 0 2]; %m -> desired space: [-4 4; -4 4; 0 4.5]
spacing = 0.5;
x_marker = flight_envelope(1,1):spacing:flight_envelope(1,2);
y_marker = flight_envelope(2,1):spacing:flight_envelope(2,2);
z_marker = flight_envelope(3,1):spacing:flight_envelope(3,2);
[X,Y,Z] = meshgrid(x_marker, y_marker, z_marker);

specs.Target = [X(:), Y(:), Z(:)];
specs.NumPoints = size(specs.Target, 1);

% Pre-compute static Data 
% section centres based on flight enveloope size and number of cameras
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

% Pre-compute data for resolution uncertainty
specs.PreComputed.adjacentSurfaces = [1 2; 2 3; 3 4; 4 1];
specs.PreComputed.du = 0.5; %pixels
specs.PreComputed.dv = 0.5; %pixels
specs.PreComputed.penaltyUncertainty = 100;
specs.PreComputed.w2 = 0.2; %weight for ellipsoid optimization

% Pre-compute data for dynamic occlusion
specs.PreComputed.minTriangAngle = 40; % degrees
specs.PreComputed.maxTriangAngle = 140; % degrees
specs.PreComputed.maxCameraRange = 16; %m effective range

% Pre-allocate camera intrinsic matrix (same for all cameras)
specs.K = [specs.Focal/specs.PixelSize, 0, specs.PrincipalPoint(1);
           0, specs.Focal/specs.PixelSize, specs.PrincipalPoint(2);
           0, 0, 1];

% Pre-compute target space in homogeneous coordinates for batch processing
specs.TargetHomogeneous = [specs.Target'; ones(1, specs.NumPoints)];

specs.PreComputed.uncertNorm = 100;  
specs.PreComputed.occlNorm = 440;    
%% Problem Definition
% Cost function:
% 1 = Resolution Uncertainty only
% 2 = Dynamic Occlusion only  
% 3 = Combined (weighted)

costFunctionType = 3; % Change this to select cost function

%Weights for combined cost function (only used if costFunctionType = 3)
%Weights need to sum to 1
specs.WeightUncertainty = 0.5; % Resolution uncertainty weight
specs.WeightOcclusion = 0.5; % Dynamic occlusion weight

% Select cost function
switch costFunctionType
    case 1
        problem.CostFunction = @resUncertainty; %moved camera chromosome to combined cost function 
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
params.MaxIt = 100;
params.nPop = 100;
params.beta = 1;
params.pC = 1;
params.gamma = 0.1;
params.mu = 0.1; %probability of mutation
params.sigma = 0.1;
params.Tournamentsize =3;

% profile on
% % Single cost function evaluation timing
% chromosome = initialPopulation(problem.VarMin, problem.VarMax, specs.SectionCentres, specs.Cams);
% 
% fprintf('\n=== Single Cost Function Breakdown ===\n');
% 
% % Time setupCameras
% tic;
% [cameras, CamCenters] = setupCameras(chromosome, specs.Cams, specs.Resolution, ...
%     specs.PixelSize, specs.Focal, specs.PrincipalPoint, specs.npix);
% t_setup = toc;
% fprintf('setupCameras: %.4f sec\n', t_setup);
% 
% % Time resUncertainty
% tic;
% uncertCost = resUncertainty(specs, cameras, CamCenters);
% t_uncert = toc;
% fprintf('resUncertainty: %.4f sec\n', t_uncert);
% 
% % Time dynamicOcclusion
% tic;
% occCost = dynamicOcclusion(specs, cameras, CamCenters);
% t_occl = toc;
% fprintf('dynamicOcclusion: %.4f sec\n', t_occl);
% 
% fprintf('\nTotal single evaluation: %.4f sec\n', t_setup + t_uncert + t_occl);
% fprintf('Estimated time for 300 pop × 100 iter: %.2f hours\n', ...
%     (t_setup + t_uncert + t_occl) * 300 * 100 / 3600);
% 
% profile off
% profile viewer
%% Run GA
tic; % start timer
out = RunGA(problem, params, specs);
elapsedTime = toc; % end timer


% %% Diagnostic: Check cost component ranges
% fprintf('\n=== Cost Component Analysis ===\n');
% numSamples = 30;
% 
% uncertVals = zeros(numSamples, 1);
% occlVals = zeros(numSamples, 1);
% 
% parfor s = 1:numSamples
%     testChrom = initialPopulation(problem.VarMin, problem.VarMax, specs.SectionCentres, specs.Cams);
%     [cams, centers] = setupCameras(testChrom, specs.Cams, specs.Resolution, ...
%         specs.PixelSize, specs.Focal, specs.PrincipalPoint, specs.npix);
%     uncertVals(s) = resUncertainty(specs, cams, centers);
%     occlVals(s) = dynamicOcclusion(specs, cams, centers);
% end
% 
% fprintf('Uncertainty - Mean: %.4f, Std: %.4f, Range: [%.4f, %.4f]\n', ...
%     mean(uncertVals), std(uncertVals), min(uncertVals), max(uncertVals));
% fprintf('Occlusion   - Mean: %.4f, Std: %.4f, Range: [%.4f, %.4f]\n', ...
%     mean(occlVals), std(occlVals), min(occlVals), max(occlVals));
% fprintf('Ratio (Occl/Uncert): %.2f\n', mean(occlVals)/mean(uncertVals));

%% Results 
currentDateTime = datetime('now');
dateTimeStr = string(currentDateTime, 'yyyyMMdd_HHmmSS');
filename = sprintf('%dCams_Run_%s.mat', specs.Cams, dateTimeStr);

% Data to save
saveData.BestSolution = out.bestsol;
saveData.BestCost = out.bestsol.Cost;
saveData.CameraConfiguration = reshape(out.bestsol.Chromosome, 6, specs.Cams)';
saveData.Specifications = specs;
saveData.GAParams = params;
saveData.ConvergenceHistory = out.bestcost;
saveData.AvgCostHistory = out.avgcost;
saveData.TopTenAvgCostHistory = out.topTenAvgCost;
saveData.ElapsedTime = elapsedTime;
saveData.Timestamp = currentDateTime;

% Determine if warm-start was used (check if pop(1) was pre-set in RunGA)
warmStartUsed = true; % Set this to true when you use warm-start

%Coverage statistics
figTitle = sprintf('Camera Coverage - %d Cameras (Cost: %.4f)', specs.Cams, saveData.BestCost);
[coverageStats] = visualizeCameraCoverage(out.bestsol.Chromosome, specs, figTitle);

coveragePlotFilename = sprintf('%dCams_Run_%s_coverage.png', specs.Cams, dateTimeStr);
saveas(gcf, coveragePlotFilename);
fprintf('Coverage plot saved to: %s\n', coveragePlotFilename);
saveData.CoverageStats = coverageStats;

% Best solution configuration 
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
fprintf('Computation Time: %.2f min\n', elapsedTime/60);

% Text summary
txtFilename = sprintf('%dCams_Run_%s.txt', specs.Cams, dateTimeStr);
fid = fopen(txtFilename, 'w');
fprintf(fid, 'Genetic Algorithm Camera Placement Results\n');
fprintf(fid, '==========================================\n\n');
fprintf(fid, 'Timestamp: %s\n', char(currentDateTime));
fprintf(fid, 'Number of Cameras: %d\n', specs.Cams);
fprintf(fid, 'Best Cost (Cost Function %d): %.6f\n', costFunctionType, saveData.BestCost);
fprintf(fid, 'Cost Function Weights: %.2f Resolution Uncertainty and %.2f Dynamic Occlusion \n', specs.WeightUncertainty , specs.WeightOcclusion);
fprintf(fid, 'Workspace Size: [%.1f %.1f; %.1f %.1f; %.1f %.1f] m\n', ...
    flight_envelope(1,1), flight_envelope(1,2), ...
    flight_envelope(2,1), flight_envelope(2,2), ...
    flight_envelope(3,1), flight_envelope(3,2));
fprintf(fid, 'Computation Time: %.2f min\n', elapsedTime/60);
fprintf(fid, 'Mutation rate: %.2f \n', params.mu);
fprintf(fid, 'Tournament Size: %.2f \n', params.Tournamentsize);
fprintf(fid, 'Warm Starting Used: %s', warmStartUsed);
fprintf(fid, 'Population Size: %d', params.nPop);
fprintf(fid, 'Number of Generations: %d', params.MaxIt);
fprintf(fid, '\n==========================================\n\n');

fprintf(fid, 'Camera Configurations:\n');
fprintf(fid, '---------------------\n');
for i = 1:specs.Cams
    fprintf(fid, '\nCamera %d:\n', i);
    fprintf(fid, '  Position (m): [%.3f, %.3f, %.3f]\n', saveData.Cameras(i).Position);
    fprintf(fid, '  Orientation (rad): [%.3f, %.3f, %.3f]\n', saveData.Cameras(i).Orientation);
    fprintf(fid, '  Orientation (deg): [%.1f, %.1f, %.1f]\n', saveData.Cameras(i).OrientationDegrees);
end

fprintf(fid, '\n==========================================\n\n');
fprintf(fid, 'Camera Coverage Statistics:\n');
fprintf(fid, '---------------------------\n');
fprintf(fid, 'Total target points: %d\n', coverageStats.numPoints);
fprintf(fid, 'Points with 0 cameras: %d (%.1f%%)\n', coverageStats.zeroCameras, coverageStats.zeroCamerasPercent);
fprintf(fid, 'Points with 1 camera: %d (%.1f%%)\n', coverageStats.oneCamera, coverageStats.oneCameraPercent);
fprintf(fid, 'Points with 2+ cameras: %d (%.1f%%)\n', coverageStats.twoPlusCameras, coverageStats.twoPlusCamerasPercent);
fprintf(fid, '\nCoverage Metrics:\n');
fprintf(fid, '  Average coverage: %.2f cameras per point\n', coverageStats.avgCoverage);
fprintf(fid, '  Maximum coverage: %d cameras per point\n', coverageStats.maxCoverage);
fprintf(fid, '  Minimum coverage: %d cameras per point\n', coverageStats.minCoverage);
fprintf(fid, '  Median coverage: %.2f cameras per point\n', coverageStats.medianCoverage);

fclose(fid);
fprintf('Summary saved to: %s\n', txtFilename);

% Convergence Plot with Subplots
figure('Name', 'GA Convergence Analysis', 'Position', [100, 100, 1200, 500]);

% Subplot 1: Best Cost Convergence
subplot(1, 2, 1);
semilogy(out.bestcost, 'LineWidth', 2, 'Color', [0 0.4470 0.7410]); 
hold on;
grid on;
xlabel('Iterations');
ylabel('Best Cost (logarithmic)');
title(sprintf('Best Cost Convergence - %d Cameras', specs.Cams));
text(0.6*params.MaxIt, max(out.bestcost)*0.5, ...
    sprintf('Final Cost: %.4f\nTime: %.1fsec', saveData.BestCost, elapsedTime), ...
    'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k');
hold off;

% Subplot 2: Average Cost Evolution
subplot(1, 2, 2);
semilogy(out.avgcost, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980], 'DisplayName', 'Population Average'); 
hold on;
semilogy(out.topTenAvgCost, 'LineWidth', 2, 'Color', [0.9290 0.6940 0.1250], 'DisplayName', 'Top 10 Average');
semilogy(out.bestcost, 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410], 'LineStyle', '--', 'DisplayName', 'Best Cost');
grid on;
xlabel('Iterations');
ylabel('Cost (logarithmic)');
title('Population Cost Evolution');
legend('Location', 'best');
hold off;

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
fprintf('Camera plot saved to: %s\n', cameraPlotFilename);

%% Append to Master Log File
masterLogFile = 'GGA_RunsLog.mat';

% Create new log entry
newLogEntry.Timestamp = currentDateTime;
newLogEntry.NumCameras = specs.Cams;
newLogEntry.BestCost = saveData.BestCost;
newLogEntry.ElapsedTime = elapsedTime;
newLogEntry.CostFunctionType = costFunctionType;
newLogEntry.WarmStart = warmStartUsed;
newLogEntry.GAParams = params;
newLogEntry.RunFilename = filename;  % Reference to full results file

% Load existing log or create new one
if isfile(masterLogFile)
    load(masterLogFile, 'runLog');
    runLog(end+1) = newLogEntry;
else
    runLog = newLogEntry;
end

% Save updated log
save(masterLogFile, 'runLog');
fprintf('Run logged to: %s (Total runs: %d)\n', masterLogFile, length(runLog));

%% Animation Plot 
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

%%
% View all logged runs
viewGALog();

% View only 7-camera runs
viewGALog('NumCameras', 7);

% View only warm-start runs
viewGALog('WarmStart', true);

% Plot comparison of 7-camera runs only
plotGARuns('NumCameras', 7);

% Plot only warm-start vs cold-start for resolution uncertainty
plotGARuns('CostFunction', 3);

% Compare 7-camera, resolution uncertainty, matching GA params
plotGARuns('NumCameras', 7, 'CostFunction', 1, 'MatchParams', true);

%%

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

