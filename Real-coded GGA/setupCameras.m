function cameras = setupCameras(cameraChromosome, numCams, resolution, pixelSize, focalLength, PrincipalPoint)
%Initialize camera objects from chromosome
    cameras = cell(numCams, 1);
    for i = 1:numCams
        chromStartIdx = (i-1)*6 + 1;
        chromEndIdx = i*6;
        camPosition = cameraChromosome(chromStartIdx:chromStartIdx+2);
        camOrientation = cameraChromosome(chromEndIdx-2:chromEndIdx);
        
        T = se3(eul2rotm(camOrientation, "XYZ"), camPosition);
        cameras{i} = CentralCamera(name="cam"+i, resolution=resolution, pixel=pixelSize, focal=focalLength, pose=T, center=PrincipalPoint);
    end
end