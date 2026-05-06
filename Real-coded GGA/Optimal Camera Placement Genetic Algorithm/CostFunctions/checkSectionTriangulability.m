function triangulable = checkSectionTriangulability(camViewVectors, frontCams, minAngle, maxAngle)
%Check if any pair of front cameras can triangulate
    triangulable = false;
    numFront = length(frontCams);
    
    for i = 1:numFront
        for j = i+1:numFront
            if checkTriangulability(camViewVectors, frontCams(i), frontCams(j), minAngle, maxAngle)
                triangulable = true;
                return;
            end
        end
    end
end