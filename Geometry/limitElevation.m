function [chrom, nChanged] = limitElevation(chrom, minElevDeg)
%LIMITELEVATION  Tilt cameras up so no optical axis points below minElevDeg.
%
%   [chrom, nChanged] = limitElevation(chrom, minElevDeg)
%
%   minElevDeg is the lowest allowed elevation of the optical axis (camera
%   +z) above horizontal [deg]. The OpenMV cameras are screwed into the
%   ball head through the back of the camera, so the 90-degree slot stops
%   them at level: they can look level or up, never down (minElevDeg = 0,
%   cfg.Mount.MinElevation).
%
%   A camera that points too low is rotated about the horizontal axis
%   perpendicular to its optical axis until it sits exactly at minElevDeg.
%   This is the smallest rotation that makes it legal: azimuth is kept and
%   roll about the optical axis is unchanged. Cameras already legal are
%   untouched. Position genes are unchanged. No-op when minElevDeg is
%   empty or <= -90.
%
%   This is a constraint on the direction R(:,3), which depends on all
%   three rotation-vector genes together, so it cannot be written as box
%   bounds (VarMin/VarMax) on the genes. It is applied after mutation,
%   like the mount regions (enforceMount).

    nChanged = 0;
    if nargin < 2 || isempty(minElevDeg) || minElevDeg <= -90
        return;
    end
    sMin = sind(minElevDeg);
    cMin = cosd(minElevDeg);
    for s = 0:6:numel(chrom) - 6
        R = genesToRotm(chrom(s+4:s+6));
        d = R(:, 3);
        if d(3) >= sMin - 1e-12, continue; end
        h = d(1:2);
        if norm(h) < 1e-9                % straight down: use the image "up" side
            h = -R(1:2, 2);
            if norm(h) < 1e-9, h = [1; 0]; end
        end
        h = h / norm(h);
        t = [cMin*h; sMin];               % same azimuth, at the limit
        ax = cross(d, t);
        ang = atan2(norm(ax), dot(d, t));
        if norm(ax) < 1e-12               % d = -t cannot occur (d below, t at/above)
            continue;
        end
        chrom(s+4:s+6) = rotmToGenes(genesToRotm(ax / norm(ax) * ang) * R);
        nChanged = nChanged + 1;
    end
end
