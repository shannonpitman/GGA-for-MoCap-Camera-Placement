function occlusionError = dynamicOcclusion(cameraChromosome, specs)
%   Based on probabilistic model from Rahimian & Kearney 2017 paper
%   Considers all possible orientations of vertical occluder
% occlusionError= Mean occlusion angle across all target points

    numCams = specs.Cams;
    resolution = specs.Resolution;
    pixelSize = specs.PixelSize;
    focalLength = specs.Focal;
    PrincipalPoint = specs.PrincipalPoint;
    TargetSpace = specs.Target;
    Npix = specs.npix;

    % Triangulability constraints (from paper)
    minTriangAngle = 40; % degrees
    maxTriangAngle = 140; % degrees
    maxCameraRange = 700; % cm effective range -> update when tested 
    
    % Initialize cameras
    cameras = setupCameras(cameraChromosome, numCams, resolution, pixelSize, focalLength, PrincipalPoint,Npix);

    numPoints = size(TargetSpace, 1);
    occlusionAngles = zeros(numPoints, 1);
    
    % Process each target point
    parfor p = 1:numPoints
        point = TargetSpace(p, :);
        
        % Find which cameras can see this point
        [visibleCams, camViewVectors] = findVisibleCameras(point, cameras, numCams, resolution, maxCameraRange);
        
        % Calculate occlusion angle for this point
        occlusionAngles(p) = calculatePointOcclusion(visibleCams, camViewVectors, minTriangAngle, maxTriangAngle);
    end
    
    occlusionError = mean(occlusionAngles)/3.6; %scale to be a point score out of 100
end





