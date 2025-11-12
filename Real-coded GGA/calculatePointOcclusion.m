function Q = calculatePointOcclusion(visibleCams, camViewVectors, minAngle, maxAngle)
% Calculate occlusion angle for a single point
% Q= Sum of angles for which target point is not visible in two triangulable views

    numVisible = length(visibleCams);
    
    %Need at least 2 cameras to triangulate
    if numVisible < 2
        Q = 360; % Maximum error
        return;
    end
    
    horizViewVectors = camViewVectors(1:2, :); % Project view vectors onto horizontal plane
    norms = vecnorm(horizViewVectors);
    validIdx = norms > 0;
    horizViewVectors(:, validIdx) = horizViewVectors(:, validIdx) ./ norms(validIdx);
    
    viewAngles = atan2d(horizViewVectors(2, :), horizViewVectors(1, :)); % Calculate angles in horizontal plane [degrees]
    viewAngles = mod(viewAngles, 360); % Convert to 0-360 range
    viewAngles = sort(viewAngles);

    Q = calculateOccludedSections(viewAngles, camViewVectors, minAngle, maxAngle); % Process each section between camera views
end
