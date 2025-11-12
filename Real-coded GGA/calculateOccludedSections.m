function Q = calculateOccludedSections(viewAngles, camViewVectors, minAngle, maxAngle)
%Calculate total angle of occluded sections
    
    numCams = length(viewAngles);
    Q = 0;
    
    %Sections between consecutive camera view angles
    for s = 1:numCams
        if s < numCams
            sectionStart = viewAngles(s);
            sectionEnd = viewAngles(s+1);
        else
            sectionStart = viewAngles(s);
            sectionEnd = viewAngles(1) + 360;
        end
        
        sectionSize = mod(sectionEnd - sectionStart, 360);
        
        %Cameras on front side of occluder for this section
        testAngle = mod(sectionStart + sectionSize/2, 360);
        frontCams = findFrontCameras(viewAngles, testAngle);
        
        %Front cameras satisfying triangulability
        triangulable = checkSectionTriangulability(camViewVectors, frontCams, minAngle, maxAngle);
        
        %If not triangulable, section angle added to error
        if ~triangulable
            Q = Q + sectionSize;
        end
    end
end
