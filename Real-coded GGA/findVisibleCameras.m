function [visibleCams, camViewVectors] = findVisibleCameras(point, cameras, camCenters, numCams, resolution, maxRange)
%Returns indices of cameras that can see the point and their view vectors
visibleCams = zeros(1, numCams);
camViewVectors = zeros(3, numCams);
visCount = 0;  % Counter for visible cameras
for i = 1:numCams
% Project point to camera image plane
   uv = cameras{i}.project(point);
   u = uv(1);
   v = uv(2);
    % Check if point is within field of view
    if (u >= 1 && u <= resolution(1) && v >= 1 && v <= resolution(2))
    % Calculate view vector from camera to point
        camCenter = cameras{i}.center();
        camCenter = camCenter(:)';
        viewVector = point - camCenter;
        distance = norm(viewVector);
        % Check if within effective range
        if distance <= maxRange && distance > 0
            visCount = visCount + 1;
            visibleCams(visCount) = i;
            camViewVectors(:, visCount) = viewVector / distance;  % 
        end
    end
    % Trim arrays to actual size
    visibleCams = visibleCams(1:visCount);
    camViewVectors = camViewVectors(:, 1:visCount);
end
