function [cameras, cameraCenters] = setupCameras(cameraChromosome, numCams, resolution, pixelSize, focalLength, PrincipalPoint, Npix)
%Initialize camera objects from chromosome
    cameras = cell(numCams, 1);
    cameraCenters = zeros(3, numCams);
    for i = 1:numCams
        chromStartIdx = (i-1)*6 + 1;
        chromEndIdx = i*6;
        camPosition = cameraChromosome(chromStartIdx:chromStartIdx+2);
        camOrientation = cameraChromosome(chromEndIdx-2:chromEndIdx);
        
        T = se3(eul2rotm(camOrientation, "XYZ"), camPosition);
        cameras{i} = CentralCamera(name="cam"+i, resolution=resolution, pixel=pixelSize, focal=focalLength, pose=T, center=PrincipalPoint, npix=Npix);
        cameraCenters(:, i) = camPosition(:);
    end
end