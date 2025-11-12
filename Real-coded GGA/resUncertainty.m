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
    
    
    adjacentSurfaces = [1 2; 2 3; 3 4; 4 1];
    du = 0.5; %pixels
    dv = 0.5; %pixels
    penaltyUncertainty = 100; %high uncertainty for points camera cannot see 
    
    %Compute Camera Transforms -> RVC3 Central Camera 
    cameras = setupCameras(cameraChromosome, numCams, resolution, pixelSize, focalLength, PrincipalPoint);    
        
    numPoints = size(TargetSpace,1);
    uncertainties = zeros(numPoints,1);
    
    parfor p =1:numPoints
        point = TargetSpace(p,:);
        uncertainties(p) = computePointUncertainty(point, cameras, numCams, resolution, adjacentSurfaces, du,dv, penaltyUncertainty, w2)
    end

    errorVolume = mean(uncertainties);
end
