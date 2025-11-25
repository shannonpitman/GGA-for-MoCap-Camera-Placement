function errorVolume = resUncertainty(specs, cameras, CamCenters)
%Computes total uncertainty due to image quantisation over a 3D target
%Space for a given chromosome (camera arrangement) 
    numCams = specs.Cams;

    %Camera Parameters
    resolution = specs.Resolution;
    TargetSpace = specs.Target;

    adjacentSurfaces = specs.PreComputed.adjacentSurfaces;
    du = specs.PreComputed.du;
    dv = specs.PreComputed.dv;
    penaltyUncertainty = specs.PreComputed.penaltyUncertainty;
    w2 = specs.PreComputed.w2;

    numPoints = specs.NumPoints;
    uncertainties = zeros(numPoints,1);
    
    parfor p =1:numPoints
        point = TargetSpace(p,:);
        uncertainties(p) = computePointUncertainty(point, cameras, CamCenters, numCams, adjacentSurfaces, du,dv, penaltyUncertainty, w2, resolution)
    end

    errorVolume = mean(uncertainties);
end
