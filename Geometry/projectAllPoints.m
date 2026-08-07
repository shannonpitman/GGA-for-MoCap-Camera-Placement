function [U, V, visMask] = projectAllPoints(cameras, TargetSpace, camCentres, resolution)
% Batched pinhole projection of every target point through every camera. 
    numCams = numel(cameras);
    N = size(TargetSpace, 1);
    Xw = TargetSpace.'; % 3 x N  (points as columns)

    U = zeros(N, numCams);
    V = zeros(N, numCams);
    visMask = false(N, numCams);

    W = resolution(1);
    H = resolution(2);

    for i = 1:numCams
        K   = cameras{i}.K;            % 3 x 3 intrinsics
        Rcw = cameras{i}.T.rotm;       % camera-to-world rotation
        c   = camCentres(:, i);        % 3 x 1 camera centre

        Xc  = Rcw.' * (Xw - c);        % 3 x N : world -> camera frame
        uvw = K * Xc;                  % 3 x N : homogeneous image coords

        depth = uvw(3, :).';              % N x 1 : projective depth
        u = (uvw(1, :) ./ uvw(3, :)).';   % N x 1
        v = (uvw(2, :) ./ uvw(3, :)).';   % N x 1

        U(:, i) = u;
        V(:, i) = v;
        % Visible only if IN FRONT of the camera (depth > 0) AND within the
        % image. The depth test matters: a point behind the camera projects
        % to a mirrored (u,v) that can still fall inside the image bounds,
        % which CentralCamera.project rejects. Omitting it over-counts
        % visibility and mis-scores those points.
        visMask(:, i) = (depth > 0 & u >= 1 & u <= W & v >= 1 & v <= H);
    end
end
