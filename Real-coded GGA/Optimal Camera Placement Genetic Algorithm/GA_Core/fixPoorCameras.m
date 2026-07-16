function y = fixPoorCameras(x, specs, coverageThreshold)
    % Quick fix for cameras that see too few points
    % coverageThreshold: minimum fraction of points a camera should see
    
    numCams = specs.Cams;
    numPoints = size(specs.Target, 1);
    minPointsRequired = ceil(coverageThreshold * numPoints);

    % Vectorised per-camera coverage, built directly from the chromosome
    % (no CentralCamera objects, no per-point loop). Replaces the old
    % setupCameras + per-point findVisibleCameras loop, which was the
    % dominant serial cost when the GA calls this thousands of times.
    cameraCoverage = cameraCoverageCounts(x, specs);

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