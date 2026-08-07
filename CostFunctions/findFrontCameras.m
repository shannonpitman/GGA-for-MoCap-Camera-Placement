function frontCams = findFrontCameras(viewAngles, testAngle)
%FINDFRONTCAMERAS Find cameras on front side of occluder
    camAngles = mod(viewAngles, 360);
    angleDiffs = mod(testAngle - camAngles, 360);
    
    % Cameras on front side have angle difference in (90, 270)
    frontCams = find(angleDiffs > 90 & angleDiffs < 270);
end