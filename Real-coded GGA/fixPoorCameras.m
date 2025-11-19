function y = fixPoorCameras(x, specs, coverageThreshold)
    % Quick fix for cameras that see too few points
    % coverageThreshold: minimum fraction of points a camera should see (e.g., 0.05 = 5%)
    
    numCams = specs.Cams;
    TargetSpace = specs.Target;
    numPoints = size(TargetSpace, 1);
    minPointsRequired = ceil(coverageThreshold * numPoints);
    
    % Setup cameras
    cameras = setupCameras(x, numCams, specs.Resolution, specs.PixelSize, specs.Focal, specs.PrincipalPoint);
    
    % Quick coverage check for each camera
    cameraCoverage = zeros(numCams, 1);
    for c = 1:numCams
        for p = 1:numPoints
            uv = cameras{c}.project(TargetSpace(p, :));
            if (uv(1) >= 1 && uv(1) <= specs.Resolution(1) && uv(2) >= 1 && uv(2) <= specs.Resolution(2))
                cameraCoverage(c) = cameraCoverage(c) + 1;
            end
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