function coverage = cameraCoverageCounts(x, specs)
%CAMERACOVERAGECOUNTS  Number of target points each camera sees.
%
%   coverage = cameraCoverageCounts(x, specs)  ->  numCams x 1 integer counts
%
%   Fully vectorised: the projection is built straight from the chromosome
%   (intrinsics K from focal/pixel size, rotation R from the Euler angles,
%   centre from the position), with NO CentralCamera objects and NO
%   per-point loop. This is a drop-in for the old
%   setupCameras + per-point findVisibleCameras coverage loop, which
%   dominated wall-clock time because fixPoorCameras runs thousands of times
%   per GA run on the (serial) client.
%
%   A (point, camera) pair counts as covered when the point is in front of
%   the camera, projects inside the image, and is within the camera's
%   effective range -- identical to findVisibleCameras. Verified against the
%   object-based path in testFixPoorCamerasEquivalence.m.

    numCams = specs.Cams;
    Xw      = specs.Target.';                 % 3 x N (points as columns)

    W  = specs.Resolution(1);
    H  = specs.Resolution(2);
    u0 = specs.PrincipalPoint(1);
    v0 = specs.PrincipalPoint(2);
    rho = specs.PixelSize;

    fN = specs.Focal;
    fW = specs.FocalWide;
    rN = specs.PreComputed.maxCameraRange;
    rW = specs.PreComputed.maxCameraRangeWide;

    coverage = zeros(numCams, 1);

    for i = 1:numCams
        s      = (i-1)*6;
        camPos = x(s+1:s+3);
        camEul = x(s+4:s+6);

        if mod(i,2) == 0            % even cameras use the wide lens (setupCameras)
            f = fW;  effRange = rW;
        else
            f = fN;  effRange = rN;
        end

        Rcw = eul2rotm(camEul, "XYZ");                 % camera-to-world (== setupCameras)
        K   = [f/rho, 0, u0; 0, f/rho, v0; 0, 0, 1];   % pinhole intrinsics
        c   = camPos(:);                               % 3 x 1 camera centre

        d   = Xw - c;                  % 3 x N : point - centre (world)
        uvw = K * (Rcw.' * d);         % 3 x N : homogeneous image coords

        depth = uvw(3, :);
        u = uvw(1, :) ./ depth;
        v = uvw(2, :) ./ depth;
        dist = sqrt(sum(d.^2, 1));     % 1 x N : range to point

        vis = depth > 0 & u >= 1 & u <= W & v >= 1 & v <= H & ...
              dist <= effRange & dist > 0;

        coverage(i) = sum(vis);
    end
end
