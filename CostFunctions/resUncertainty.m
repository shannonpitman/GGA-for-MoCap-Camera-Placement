function [errorVolume, uncertainties] = resUncertainty(specs, cameras, CamCenters)
% Second output: per-point uncertainty [m], for cost-field plots.
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
    % loop. Visibility is the shared FOV + range test used by the occlusion
    % term and the coverage statistics. Returns U, V, visMask (numPoints x numCams).
    [visMask, ~, U, V] = projectVisibilityOcclusion(cameras, TargetSpace, CamCenters, ...
        resolution, specs.PreComputed.maxCameraRange, specs.PreComputed.maxCameraRangeWide, ...
        specs.FocalWide);

    parfor p =1:numPoints
        uncertainties(p) = computePointUncertainty(TargetSpace(p,:), cameras, CamCenters, ...
            numCams, adjacentSurfaces, du,dv, penaltyUncertainty, w2, resolution, ...
            U(p,:), V(p,:), visMask(p,:));
    end

    errorVolume = mean(uncertainties);
end
