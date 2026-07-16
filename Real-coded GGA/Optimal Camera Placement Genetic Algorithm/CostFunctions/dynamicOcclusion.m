function occlusionError = dynamicOcclusion(specs, cameras, CamCenters)
%   Based on probabilistic model from Rahimian & Kearney 2017 paper
%   Considers all possible orientations of vertical occluder
% occlusionError= Mean occlusion angle across all target points

    numCams = specs.Cams;
    resolution = specs.Resolution;
    TargetSpace = specs.Target;

    % Triangulability constraints (from paper)
    minTriangAngle = specs.PreComputed.minTriangAngle; % degrees
    maxTriangAngle = specs.PreComputed.maxTriangAngle; % degrees
    maxCameraRange = specs.PreComputed.maxCameraRange; %m -> 16m effective range (OptiTrack)-> update when tested 
    maxCameraRangeWide = specs.PreComputed.maxCameraRangeWide;
    focalWide = specs.FocalWide;
    numPoints = size(TargetSpace, 1);
    occlusionAngles = zeros(numPoints, 1);

    % Batched visibility + view vectors for every point/camera, computed
    % once here instead of numPoints*numCams per-point project() + range
    % calls inside findVisibleCameras.
    [visMask, viewUnit] = projectVisibilityOcclusion(cameras, TargetSpace, ...
        CamCenters, resolution, maxCameraRange, maxCameraRangeWide, focalWide);

    % Process each target point
    parfor p = 1:numPoints
        visRow  = visMask(p, :);
        viewMat = reshape(viewUnit(p, :, :), numCams, 3).';   % 3 x numCams

        % Find which cameras can see this point (from precomputed data)
        [visibleCams, camViewVectors] = findVisibleCameras(TargetSpace(p, :), cameras, ...
            CamCenters, numCams, resolution, maxCameraRange, maxCameraRangeWide, focalWide, ...
            visRow, viewMat);

        % Calculate occlusion angle for this point
        occlusionAngles(p) = calculatePointOcclusion(visibleCams, camViewVectors, minTriangAngle, maxTriangAngle);
    end
    
    occlusionError = mean(occlusionAngles); %scale to be a point score out of 100
end





