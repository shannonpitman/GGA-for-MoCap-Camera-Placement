function visualizeCameraCoverage(cameraChromosome, specs, figTitle)
%Visualizes camera coverage of target points with color coding:
% Green: No cameras see the point
% Magenta: Only one camera sees the point  
% Red: Two or more cameras see the point

    numCams = specs.Cams;
    
    %Camera Parameters
    resolution = specs.Resolution;
    pixelSize = specs.PixelSize;
    focalLength = specs.Focal;
    PrincipalPoint = specs.PrincipalPoint;
    TargetSpace = specs.Target;
    
    cameras = cell(numCams,1);
    
    %Compute Camera Transforms
    for i = 1:numCams
        chromStartIdx = (i-1)*6+1;
        chromEndIdx = i*6;
        camPositions = cameraChromosome(chromStartIdx: chromStartIdx+2);
        camOrientations = cameraChromosome(chromEndIdx-2: chromEndIdx);
    
        T = se3(eul2rotm(camOrientations, "XYZ"), camPositions);
        cameras{i} = CentralCamera(name="cam"+i, resolution=resolution, pixel=pixelSize, focal=focalLength, pose=T, center=PrincipalPoint); 
    end
    
    numPoints = size(TargetSpace,1);
    cameraCoverage = zeros(numPoints,1);
    
    %Amount of cams seeing each point
    for p = 1:numPoints
        point = TargetSpace(p,:);
        visibleCount = 0;
        
        for i = 1:numCams
            uv = cameras{i}.project(point);
            u = uv(1);
            v = uv(2);
            
            %In camera's field of view
            if (u >= 1 && u <= resolution(1) && v >= 1 && v <= resolution(2))
                visibleCount = visibleCount + 1;
            end
        end
        
        cameraCoverage(p) = visibleCount;
    end
    
    % Colour map based on coverage
    colors = zeros(numPoints, 3);
    parfor p = 1:numPoints
        if cameraCoverage(p) == 0
            colors(p,:) = [0, 1, 0]; % Green: no cameras
        elseif cameraCoverage(p) == 1
            colors(p,:) = [1, 0, 1]; % Magenta: one camera
        else
            colors(p,:) = [1, 0, 0]; % Red: two or more cameras
        end
    end
    
    %Colour coded points
    figure('Name', figTitle, 'Position', [100, 100, 800, 600]);
    scatter3(TargetSpace(:,1), TargetSpace(:,2), TargetSpace(:,3), 30, colors, 'filled');
    
    % Add cameras to the plot
    % hold on;
    % for i = 1:numCams
    %     hold on;
    %     cameras{i}.plot_camera('label', 'scale', 0.3);
    % end
    
    axis equal;
    grid on;
    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('Z (m)');
    title(figTitle);
    view(45, 30);
    
    % Add legend
    % h1 = scatter3(nan, nan, nan, 50, [0, 1, 0], 'filled');
    % h2 = scatter3(nan, nan, nan, 50, [0, 1, 1], 'filled');
    % h3 = scatter3(nan, nan, nan, 50, [1, 0, 0], 'filled');
    % legend([h1, h2, h3], {'No cameras', '1 camera', '2+ cameras'}, ...
    %     'Location', 'best');
    
    % Print statistics
    fprintf('\n Camera Coverage Statistics \n');
    fprintf('Points with 0 cameras: %d (%.1f%%)\n', ...
        sum(cameraCoverage == 0), 100*sum(cameraCoverage == 0)/numPoints);
    fprintf('Points with 1 camera: %d (%.1f%%)\n', ...
        sum(cameraCoverage == 1), 100*sum(cameraCoverage == 1)/numPoints);
    fprintf('Points with 2+ cameras: %d (%.1f%%)\n', ...
        sum(cameraCoverage >= 2), 100*sum(cameraCoverage >= 2)/numPoints);
    
    % Calculate average coverage
    avgCoverage = mean(cameraCoverage);
    fprintf('Average camera coverage: %.2f cameras per point\n', avgCoverage);
    
    hold off;
end
