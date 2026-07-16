function [U, V, visMask] = projectAllPoints(cameras, TargetSpace, camCentres, resolution)
%PROJECTALLPOINTS  Batched pinhole projection of every target point through
%   every camera. Vectorised replacement for the per-point
%   cameras{i}.project(point) calls that previously ran inside
%   computePointUncertainty (and can be reused by findVisibleCameras).
%
%   [U, V, visMask] = projectAllPoints(cameras, TargetSpace, camCentres, resolution)
%
%   Inputs:
%     cameras     - numCams x 1 (or 1 x numCams) cell of CentralCamera objects
%     TargetSpace - N x 3 world points (specs.Target)
%     camCentres  - 3 x numCams camera centres (CamCenters from setupCameras)
%     resolution  - [W H] image resolution
%
%   Outputs (all N x numCams):
%     U, V     - projected pixel coordinates for each (point, camera)
%     visMask  - true where the projection lands within [1,W] x [1,H]
%
%   Model: p ~ K * R_cw.' * (Xw - c). This is the exact forward inverse of
%   the back-projection used in quantToWorld
%   (Xw = R_cw*(K\[u;v;1]) + c), so U/V match CentralCamera.project to
%   machine precision. Verified by testProjectAllPointsEquivalence.m.
%
%   The whole point loop is replaced by one matrix multiply per camera, so
%   projection work drops from numPoints*numCams object-method calls to
%   numCams BLAS multiplies.

    numCams = numel(cameras);
    N       = size(TargetSpace, 1);
    Xw      = TargetSpace.';           % 3 x N  (points as columns)

    U       = zeros(N, numCams);
    V       = zeros(N, numCams);
    visMask = false(N, numCams);

    W = resolution(1);
    H = resolution(2);

    for i = 1:numCams
        K   = cameras{i}.K;            % 3 x 3 intrinsics
        Rcw = cameras{i}.T.rotm;       % camera-to-world rotation
        c   = camCentres(:, i);        % 3 x 1 camera centre

        Xc  = Rcw.' * (Xw - c);        % 3 x N : world -> camera frame
        uvw = K * Xc;                  % 3 x N : homogeneous image coords

        u = (uvw(1, :) ./ uvw(3, :)).';   % N x 1
        v = (uvw(2, :) ./ uvw(3, :)).';   % N x 1

        U(:, i) = u;
        V(:, i) = v;
        visMask(:, i) = (u >= 1 & u <= W & v >= 1 & v <= H);
    end
end
