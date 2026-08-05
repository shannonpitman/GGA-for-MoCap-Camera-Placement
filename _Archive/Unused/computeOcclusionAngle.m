function Q = calculatePointOcclusion(visibleCams, camViewVectors, minAngle, maxAngle)
%CALCULATEPOINTOCCLUSION Calculate occlusion angle for a single point
%   Q: Sum of angles for which target point is not visible in two triangulable views

    numVisible = length(visibleCams);
    
    % Need at least 2 cameras to triangulate
    if numVisible < 2
        Q = 360; % Maximum error
        return;
    end
    
    % Project view vectors onto horizontal plane
    horizViewVectors = camViewVectors(1:2, :);
    
    % Normalize horizontal components
    norms = vecnorm(horizViewVectors);
    validIdx = norms > 0;
    horizViewVectors(:, validIdx) = horizViewVectors(:, validIdx) ./ norms(validIdx);
    
    % Calculate angles in horizontal plane
    viewAngles = atan2d(horizViewVectors(2, :), horizViewVectors(1, :));
    viewAngles = mod(viewAngles, 360); % Convert to 0-360 range
    viewAngles = sort(viewAngles);
    
    % Process each section between camera views
    Q = calculateOccludedSections(viewAngles, camViewVectors, minAngle, maxAngle);
end
