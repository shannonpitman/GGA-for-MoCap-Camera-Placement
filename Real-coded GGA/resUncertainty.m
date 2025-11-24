function errorVolume = resUncertainty(cameraChromosome, specs)
%Computes total uncertainty due to image quantisation over a 3D target
%Space for a given chromosome (camera arrangement) 
    numCams = specs.Cams;

    %Camera Parameters
    resolution = specs.Resolution;
    pixelSize = specs.PixelSize;
    focalLength = specs.Focal;
    PrincipalPoint = specs.PrincipalPoint;
    TargetSpace = specs.Target;
    npix = specs.npix;
    
    adjacentSurfaces = specs.PreComputed.adjacentSurfaces;
    du = specs.PreComputed.du;
    dv = specs.PreComputed.dv;
    penaltyUncertainty = specs.PreComputed.penaltyUncertainty;
    w2 = specs.PreComputed.w2;
    %Compute Camera Transforms -> RVC3 Central Camera 
    [cameras, cameraCenters] = setupCameras(cameraChromosome, numCams, resolution, pixelSize, focalLength, PrincipalPoint, npix);    
        
    numPoints = specs.NumPoints;
    uncertainties = zeros(numPoints,1);
    
    parfor p =1:numPoints
        point = TargetSpace(p,:);
        uncertainties(p) = computePointUncertainty(point, cameras, cameraCenters, numCams, adjacentSurfaces, du,dv, penaltyUncertainty, w2)
    end

    errorVolume = mean(uncertainties);
end
