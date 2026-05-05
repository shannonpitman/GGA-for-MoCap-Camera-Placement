function y = fixPoorCameras(x, specs, coverageThreshold)
    % Quick fix for cameras that see too few points
    % coverageThreshold: minimum fraction of points a camera should see
    
    numCams = specs.Cams;
    TargetSpace = specs.Target;
    numPoints = size(TargetSpace, 1);
    minPointsRequired = ceil(coverageThreshold * numPoints);
    
    % Setup cameras
    resolution = specs.Resolution;
    focalLength = specs.Focal;
    focalLengthWide = specs.FocalWide;
    PrincipalPoint = specs.PrincipalPoint;
    maxRange = specs.PreComputed.maxCameraRange;
    maxRangeWide = specs.PreComputed.maxCameraRangeWide;
 
    [cameras, camCenters] = setupCameras(x, numCams, resolution, focalLength, focalLengthWide, PrincipalPoint);
    
    % Coverage check for each camera (with range check via findVisibleCameras)
    cameraCoverage = zeros(numCams, 1);
    for p = 1:numPoints
        point = TargetSpace(p, :);
        [visibleCams, ~] = findVisibleCameras(point, cameras, camCenters, numCams, resolution, maxRange, maxRangeWide, focalLengthWide);
        for k = 1:length(visibleCams)
            cameraCoverage(visibleCams(k)) = cameraCoverage(visibleCams(k)) + 1;
        end
    end
    
    % Fix cameras that see too few points
    y = x;
    for c = 1:numCams
        if cameraCoverage(c) < minPointsRequired
            chromStart = (c-1)*6 + 1;
            chromEnd = c*6;
            
            % Get current camera position
            camPos = x(chromStart:chromStart+2);
            
            % Find nearest section center
            distances = vecnorm(specs.SectionCentres - camPos, 2, 2);
            [~, closestIdx] = min(distances);
            targetPoint = specs.SectionCentres(closestIdx, :);
            
            % Reorient camera toward the target
            directionVec = targetPoint - camPos;
            directionUnit = directionVec / norm(directionVec);
            
            Z_axis = [0, 0, 1];
            rotAng = acos(dot(Z_axis, directionUnit));
            rotAxis = cross(Z_axis, directionUnit);
            
            if norm(rotAxis) > 0
                rotAxis = rotAxis / norm(rotAxis);
                q = quaternion([cos(rotAng/2), rotAxis * sin(rotAng/2)]);
                q = normalize(q);
                q2Eul = euler(q, 'XYZ', 'point');
                
                % Update only the orientation (keep position)
                y(chromEnd-2:chromEnd) = q2Eul;
            end
        end
    end
end