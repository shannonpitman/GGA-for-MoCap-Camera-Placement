function [pos, regionIdx] = projectToMount(pos, regions)
%PROJECTTOMOUNT  Nearest mountable position to a camera position.
%
%   [pos, regionIdx] = projectToMount(pos, regions) moves the 1x3 position to
%   the closest point in any region from mountRegions. Orientation is not
%   touched. With no regions the position is returned unchanged.

    regionIdx = 0;
    if isempty(regions)
        return;
    end

    best = inf;
    bestPos = pos;
    for k = 1:numel(regions)
        cand = nearestInRegion(pos, regions(k));
        d = sum((cand - pos).^2);
        if d < best
            best = d;
            bestPos = cand;
            regionIdx = k;
        end
    end
    pos = bestPos;
end

function q = nearestInRegion(p, r)
    q = min(max(p, r.Lo), r.Hi);            % clamp into the bounding box
    if strcmp(r.Type, 'ring')
        inside = q(1) > r.InnerLo(1) && q(1) < r.InnerHi(1) && ...
                 q(2) > r.InnerLo(2) && q(2) < r.InnerHi(2);
        if inside
            % Push out through the nearest edge of the excluded footprint
            gaps = [q(1) - r.InnerLo(1), r.InnerHi(1) - q(1), ...
                    q(2) - r.InnerLo(2), r.InnerHi(2) - q(2)];
            [~, e] = min(gaps);
            switch e
                case 1, q(1) = r.InnerLo(1);
                case 2, q(1) = r.InnerHi(1);
                case 3, q(2) = r.InnerLo(2);
                case 4, q(2) = r.InnerHi(2);
            end
        end
    end
end
