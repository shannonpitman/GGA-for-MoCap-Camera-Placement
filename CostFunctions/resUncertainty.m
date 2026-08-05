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

    % Batched projection of every point through every camera, done once
    % here instead of numPoints*numCams per-point project() calls inside the
    % loop. Returns U, V, visMask (all numPoints x numCams).
    [U, V, visMask] = projectAllPoints(cameras, TargetSpace, CamCenters, resolution);

    parfor p =1:numPoints
        uncertainties(p) = computePointUncertainty(TargetSpace(p,:), cameras, CamCenters, ...
            numCams, adjacentSurfaces, du,dv, penaltyUncertainty, w2, resolution, ...
            U(p,:), V(p,:), visMask(p,:));
    end

    errorVolume = mean(uncertainties);
end
