function [visibleCams, camViewVectors] = findVisibleCameras(point, cameras, numCams, resolution, maxRange)
%Returns indices of cameras that can see the point and their view vectors

    visibleCams = [];
    camViewVectors = zeros(3, 0);
    
    point = point(:)';
    
    for i = 1:numCams
        % Project point to camera image plane
        uv = cameras{i}.project(point);
        u = uv(1);
        v = uv(2);
        
        % Check if point is within field of view
        if (u >= 1 && u <= resolution(1) && v >= 1 && v <= resolution(2))
            % Get camera center
            camCenter = cameras{i}.center();
            camCenter = camCenter(:)'; % Ensure row vector
            
            % Calculate view vector from camera to point
            viewVector = point - camCenter;
            distance = norm(viewVector);
            
            % Check if within effective range
            if distance <= maxRange && distance > 0
                visibleCams(end+1) = i;
                normalizedVector = viewVector / distance; %just direction of camera to point 
                camViewVectors = [camViewVectors, normalizedVector(:)];
            end
        end
    end
end