function pos = sampleMount(regions)
%SAMPLEMOUNT  Uniform random camera position in a random mountable region.
%   Each region is equally likely (as the old one-of-five-faces rule), then
%   the position is uniform within it. Tripod-ring samples inside the
%   capture footprint are redrawn.

    r = regions(randi(numel(regions)));
    while true
        pos = r.Lo + rand(1, 3) .* (r.Hi - r.Lo);
        if ~strcmp(r.Type, 'ring') || ...
                ~(pos(1) > r.InnerLo(1) && pos(1) < r.InnerHi(1) && ...
                  pos(2) > r.InnerLo(2) && pos(2) < r.InnerHi(2))
            return;
        end
    end
end
