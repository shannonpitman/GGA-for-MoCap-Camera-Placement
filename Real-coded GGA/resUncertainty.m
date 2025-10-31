function errorVolume = resUncertainty(cameraChromosome, specs)
%Computes total uncertainty due to image quantisation over a 3D target
%Space for a given chromosome (camera arrangement) 
    numCams = specs.Cams;
    w2 = 0.2; %weight for second energy function
    
    %Camera Parameters
    resolution = specs.Resolution;
    pixelSize = specs.PixelSize;
    focalLength = specs.Focal;
    PrincipalPoint = specs.PrincipalPoint;
    TargetSpace = specs.Target;
    
    cameras = cell(numCams,1);
    
    adjacentSurfaces = [1 2; 2 3; 3 4; 4 1];
    du = 0.5; %pixels
    dv = 0.5; %pixels
    penaltyUncertainty = 100; %high uncertainty for points camera cannot see 
    
    %Compute Camera Transforms
    for i = 1:numCams
        chromStartIdx = (i-1)*6+1;
        chromEndIdx = i*6;
        camPositions = cameraChromosome(chromStartIdx: chromStartIdx+2);
        camOrientations = cameraChromosome(chromEndIdx-2: chromEndIdx);
    
        T = se3(eul2rotm(camOrientations, "XYZ"), camPositions); %camera to world cTw
        cameras{i} = CentralCamera(name="cam"+i,resolution= resolution, pixel= pixelSize, focal= focalLength, pose=T, center = PrincipalPoint); 
    end    
        
    numPoints = size(TargetSpace,1);
    uncertainties = zeros(numPoints,1);
    
    parfor p =1:numPoints
        point = TargetSpace(p,:);
        uncertainties(p) = computePointUncertainty(point, cameras, numCams, resolution, adjacentSurfaces, du,dv, penaltyUncertainty, w2)
    end
    
    errorVolume = mean(uncertainties);
end
