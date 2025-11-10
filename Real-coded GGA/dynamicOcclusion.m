function occlusionError = dynamicOcclusion(cameraChromosome, specs)
%Computes dynamic occlusion error based on probabilistic model from 
%Rahimian & Kearney 2017 paper. Considers all possible orientations of
%vertical occluder and calculates angles where point is not visible in 
%at least two triangulable views

    numCams = specs.Cams;
    
    %Camera Parameters
    resolution = specs.Resolution;
    pixelSize = specs.PixelSize;
    focalLength = specs.Focal;
    PrincipalPoint = specs.PrincipalPoint;
    TargetSpace = specs.Target;
    
    cameras = cell(numCams,1);
    
    % Triangulability constraint angles (from paper)
    minTriangAngle = 40; % degrees
    maxTriangAngle = 140; % degrees
    
    %Compute Camera Transforms
    for i = 1:numCams
        chromStartIdx = (i-1)*6+1;
        chromEndIdx = i*6;
        camPositions = cameraChromosome(chromStartIdx: chromStartIdx+2);
        camOrientations = cameraChromosome(chromEndIdx-2: chromEndIdx);
    
        T = se3(eul2rotm(camOrientations, "XYZ"), camPositions);
        cameras{i} = CentralCamera(name="cam"+i, resolution=resolution, ...
            pixel=pixelSize, focal=focalLength, pose=T, center=PrincipalPoint); 
    end
    
    numPoints = size(TargetSpace,1);
    occlusionAngles = zeros(numPoints,1);
    
    parfor p = 1:numPoints
        point = TargetSpace(p,:);
        occlusionAngles(p) = computeOcclusionAngle(point, cameras, numCams, ...
            resolution, minTriangAngle, maxTriangAngle);
    end
    
    occlusionError = mean(occlusionAngles);
end

function Q = computeOcclusionAngle(point, cameras, numCams, resolution, minAngle, maxAngle)
    % Q: Sum of angles for which target point is not visible in two triangulable views
    
    % Find cameras that can see the point
    visibleCams = [];
    camViewVectors = [];
    
    for i = 1:numCams
        uv = cameras{i}.project(point);
        u = uv(1);
        v = uv(2);
        
        % Check if point is in camera's field of view and range
        if (u >= 1 && u <= resolution(1) && v >= 1 && v <= resolution(2))
            camCenter = cameras{i}.center().';
            viewVector = point - camCenter;
            
            % Check distance (effective range ~365cm from paper)
            if norm(viewVector) <= 365
                visibleCams(end+1) = i;
                camViewVectors(:,end+1) = viewVector/norm(viewVector);
            end
        end
    end
    
    numVisible = length(visibleCams);
    
    if numVisible < 2
        Q = 360; % Maximum error if less than 2 cameras see the point
        return;
    end
    
    % Project view vectors onto horizontal plane for occlusion analysis
    % (vertical occluder rotates around vertical axis)
    horizViewVectors = camViewVectors(1:2,:); % Only x,y components
    horizViewVectors = horizViewVectors ./ vecnorm(horizViewVectors);
    
    % Calculate angles of each camera view in horizontal plane
    viewAngles = atan2d(horizViewVectors(2,:), horizViewVectors(1,:));
    viewAngles = mod(viewAngles, 360); % Convert to 0-360 range
    viewAngles = sort(viewAngles);
    
    % Calculate sections between view vectors
    numSections = numVisible;
    sectionRanges = zeros(numSections, 2);
    
    for i = 1:numSections
        if i < numSections
            sectionRanges(i,:) = [viewAngles(i), viewAngles(i+1)];
        else
            sectionRanges(i,:) = [viewAngles(i), viewAngles(1)+360];
        end
    end
    
    Q = 0;
    
    % For each section, check if triangulability constraint is satisfied
    for s = 1:numSections
        sectionSize = mod(sectionRanges(s,2) - sectionRanges(s,1), 360);
        
        % Find which cameras are on front side for this section
        % (using middle angle of section as test angle)
        testAngle = mod(sectionRanges(s,1) + sectionSize/2, 360);
        
        frontCams = [];
        for c = 1:numVisible
            camAngle = mod(viewAngles(c), 360);
            % Check if camera is on front side
            angleDiff = mod(testAngle - camAngle, 360);
            if angleDiff > 90 && angleDiff < 270
                frontCams(end+1) = c;
            end
        end
        
        % Check if any pair of front cameras satisfies triangulability
        triangulable = false;
        for i = 1:length(frontCams)
            for j = i+1:length(frontCams)
                % Calculate 3D convergence angle between cameras
                vec1 = camViewVectors(:, frontCams(i));
                vec2 = camViewVectors(:, frontCams(j));
                convergenceAngle = acosd(dot(vec1, vec2));
                
                if convergenceAngle >= minAngle && convergenceAngle <= maxAngle
                    triangulable = true;
                    break;
                end
            end
            if triangulable
                break;
            end
        end
        
        % If not triangulable, add this section's angle to error
        if ~triangulable
            Q = Q + sectionSize;
        end
    end
end