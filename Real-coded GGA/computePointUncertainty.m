function uncertainty = computePointUncertainty(point, cameras, numCams, resolution, adjacentSurfaces, du,dv, penaltyUncertainty, w2)
    isVisible = false(numCams,1); %logical vector to check whether a camera can see point
    planes = cell(1, numCams);

    for i = 1:numCams
        uv = cameras{i}.project(point); %uv projected coordinates 
        u = uv(1);
        v = uv(2);

        if (u >= 1 && u <= resolution(1) && v >= 1 && v <= resolution(2))
            cameraCentres = cameras{i}.center().'; %world location of Camera center
            worldPoints = quantToWorld(cameras{i}, u,v, du, dv, cameraCentres);
            planes{i} = buildPyramidSurf(cameraCentres, worldPoints, adjacentSurfaces);
            isVisible(i) = true;
        end
    end

    visibleIdx = find(isVisible);
    numVisible = length(visibleIdx);


    if numVisible == 0 %if there are less than 2 cameras then triangulation is not possible 
        uncertainty = penaltyUncertainty; %point is not visible in any cameras FOV
    elseif numVisible== 1 %if at least one camera sees the point
        uncertainty = 0.5*penaltyUncertainty;
    else
        vertices = calcVertices(numVisible, visibleIdx, adjacentSurfaces, planes);
        unique_vertices = unique(round(vertices,6),'rows', 'stable');

        V_ = unique_vertices-point; %centre the points at the mean 
        C0v_vert = cov(V_); %covariance matrix
        [A,~] = eig(C0v_vert); %[eigenvectors, eigenvalues in a diagonal matrix]
        
        
        transformed_v = (A*V_.').';%rotated points 

        %Initial guess
        axes0 = max(abs(transformed_v), [], 1).' ; %1x3 [a0,b0,c0], max magnitude of the x,y,z components of v 

        fun = @(axes) optimiseEllipsoid(axes, transformed_v.', axes0, w2);
        options = optimset('Display', 'off', 'TolX', 1e-6, 'TolFun', 1e-6, 'MaxIter', 1000, 'MaxFunEvals', 5000);
        [axes_opt, ~] = fminsearch(fun, axes0, options);

        %optimised axes [m]
        a = abs(axes_opt(1));
        b = abs(axes_opt(2));
        c = abs(axes_opt(3));

        uncertainty = 4/3*pi*a*b*c;%m^3
        
    end
end
