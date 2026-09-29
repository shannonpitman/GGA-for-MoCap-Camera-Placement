function y = fixPoorCameras(x, specs, coverageThreshold)
    % Quick fix for cameras that see too few points
    % coverageThreshold: minimum fraction of points a camera should see
    
    numCams = specs.Cams;
    numPoints = size(specs.Target, 1);
    minPointsRequired = ceil(coverageThreshold * numPoints);
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
            
            % Update only the orientation (keep position): axis-angle aim
            % stored directly as the rotation-vector genes
            y(chromEnd-2:chromEnd) = aimRotvec(directionUnit);
        end
    end
end