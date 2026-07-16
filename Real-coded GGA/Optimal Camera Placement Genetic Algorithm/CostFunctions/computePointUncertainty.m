function uncertainty = computePointUncertainty(point, cameras, cameraCentres, numCams, adjacentSurfaces, du,dv, penaltyUncertainty, w2, Resolution, uRow, vRow, visRow)
%   Optional trailing args uRow, vRow, visRow (each 1 x numCams) let a
%   caller supply the projection of THIS point precomputed in bulk by
%   projectAllPoints. When omitted the point is projected per-camera here,
%   so existing callers (e.g. plotCostField_GAvsOptiTrack) are unaffected.
    planes = cell(1, numCams);

    if nargin < 13 || isempty(visRow)
        % Original path: project this point through each camera locally.
        uAll = zeros(numCams,1);
        vAll = zeros(numCams,1);
        for i = 1:numCams
            uv = cameras{i}.project(point); %uv projected coordinates
            uAll(i) = uv(1);
            vAll(i) = uv(2);
        end
        visIdx = (uAll >= 1 & uAll <= Resolution(1) & vAll >= 1 & vAll <= Resolution(2));
    else
        % Precomputed (batched) projection supplied by projectAllPoints.
        uAll   = uRow(:);
        vAll   = vRow(:);
        visIdx = logical(visRow(:));
    end

    % Build the quantisation pyramid only for cameras that see the point.
    visibleIdx = find(visIdx);
    for k = 1:numel(visibleIdx)
        i = visibleIdx(k);
        worldPoints = quantToWorld(cameras{i}, uAll(i), vAll(i), du, dv, cameraCentres(:,i));
        planes{i} = buildPyramidSurf(cameraCentres(:,i), worldPoints, adjacentSurfaces);
    end

    numVisible = numel(visibleIdx);
    
    %early exit if triangulation is not possible
    if numVisible == 0 
        uncertainty = penaltyUncertainty; %point is not visible in any cameras FOV
        return;
    elseif numVisible== 1 %if at least one camera sees the point
        uncertainty = 0.5*penaltyUncertainty;
        return;
    end

    vertices = calcVertices(numVisible, visibleIdx, planes);
    if isempty(vertices)
        uncertainty = 0.5*penaltyUncertainty;
        return;
    end
    unique_vertices = unique(round(vertices,6),'rows', 'stable');
    if size(unique_vertices, 1) < 4
        %Fewer than 4 unique vertices means the polytope is degenerate in 3D (point/line/plane)
        uncertainty = 0.5*penaltyUncertainty;
        return;
    end

    V_centered = unique_vertices-point; %centre the data at the point (zero mean) -> standardize the data
    C0v_vert = cov(V_centered); %covariance matrix

    Eigs = eig(C0v_vert); %[eigenvectors, eigenvalues in a diagonal matrix]
    uncertainty = sum(sqrt(abs(Eigs))); % uncertainty based on the trace of the cov matrix
    if ~isfinite(uncertainty)
        uncertainty = 0.5*penaltyUncertainty;
    end

    % Old code that found error volume based on Wu, Sharma and Huang
    % formulation
    % [A,~] = eig(C0v_vert); %[eigenvectors, eigenvalues in a diagonal matrix] -> eigenvectors: new axes which our data lies on, eigenvalues = variances along these new axes
    % 
    % 
    % transformed_v = (A*V_centered.').';%rotated points 
    % 
    % 
    % 
    % %Initial guess
    % axes0 = max(abs(transformed_v), [], 1).' ; %1x3 [a0,b0,c0], max magnitude of the x,y,z components of v 
    % 
    % fun = @(axes) optimiseEllipsoid(axes, transformed_v.', axes0, w2);
    % options = optimset('Display', 'off', 'TolX', 1e-6, 'TolFun', 1e-6, 'MaxIter', 1000, 'MaxFunEvals', 5000);
    % [axes_opt, ~] = fminsearch(fun, axes0, options);
    % 
    % %optimised axes [m]
    % a = abs(axes_opt(1));
    % b = abs(axes_opt(2));
    % c = abs(axes_opt(3));
    % 
    % uncertainty = 4/3*pi*a*b*c;%m^3
end

