function [cameras, cameraCenters] = setupCameras(cameraChromosome, numCams, resolution, focalLength, focalLengthWide, PrincipalPoint)
%Initialize camera objects from chromosome
    cameras = cell(numCams, 1);
    cameraCenters = zeros(3, numCams);
    for i = 1:numCams
        chromStartIdx = (i-1)*6 + 1;
        chromEndIdx = i*6;
        camPosition = cameraChromosome(chromStartIdx:chromStartIdx+2);
        camOrientation = cameraChromosome(chromEndIdx-2:chromEndIdx);
        
        T = se3(eul2rotm(camOrientation, "XYZ"), camPosition);
        if mod(numCams,2) == 1 
            f = focalLengthWide;
        else
            f = focalLength;
        end
        cameras{i} = CentralCamera(name="cam"+i, resolution=resolution, focal=f, pose=T, center=PrincipalPoint);
        cameraCenters(:, i) = camPosition(:);
    end
end