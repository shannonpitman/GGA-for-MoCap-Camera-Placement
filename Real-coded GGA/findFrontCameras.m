function frontCams = findFrontCameras(viewAngles, testAngle)
%FINDFRONTCAMERAS Find cameras on front side of occluder
    frontCams = [];
    numCams = length(viewAngles);
    
    for c = 1:numCams
        camAngle = mod(viewAngles(c), 360);
        angleDiff = mod(testAngle - camAngle, 360);
        
        % Camera is on front side if angle difference is in (90, 270)
        if angleDiff > 90 && angleDiff < 270
            frontCams(end+1) = c;
        end
    end
end