function regions = mountRegions(cfg)
%MOUNTREGIONS  Mountable camera-position regions for a runConfig.
%
%   regions = mountRegions(cfg) returns a struct array with fields
%     Name   - label
%     Type   - 'plane' (wall/ceiling) or 'ring' (floor-standing tripods)
%     Lo, Hi - 1x3 lower/upper corner of the region's bounding box. For a
%              plane Lo and Hi are equal on the fixed axis.
%     InnerLo, InnerHi - 'ring' only: x-y rectangle that is excluded
%              (the capture footprint), so tripods never stand inside it.
%
%   The room is the position part of the camera search box
%   (cfg.CamLowerBounds/CamUpperBounds, genes 1-3); the capture footprint is
%   the x-y extent of cfg.Volume. cfg.Mount.Model selects the regions:
%     'walls_ceiling_tripod'  4 walls + ceiling + tripod ring (lab default)
%     'walls_ceiling'         4 walls + ceiling
%     'tripod'                tripod ring only
%     'tripod_cabinets'       tripods in the floor band on the two x sides
%                             (TripodHeight), tripods on floor or tables at
%                             the two y ends (Mount.EndHeight) and tripods on
%                             the cabinet tops behind the x-side band
%                             (Mount.Cabinet); see below
%     'box'                   no restriction beyond the search box (legacy)
%   cfg.Mount.TripodHeight = [zMin zMax] sets the tripod head height range.
%   cfg.Mount.TripodFootRadius = r (default 0) shrinks the 'tripod' and
%   'tripod_cabinets' regions so a tripod centre stays r from the room
%   walls, the cabinet fronts and the cabinet-top edges (its legs need the
%   room). The flight-space edge is not inset: feet may stand just inside
%   it, below the lowest target.
%   'tripod_cabinets' also reads
%     cfg.Mount.EndHeight  [zMin zMax] head height at the y ends (floor or table)
%     cfg.Mount.Cabinet    struct: XDepth [xIn xOut] (|x| of the cabinet top,
%                          both x sides), YSpan [y0 y1], Height (cabinet top
%                          [m]); heads at Height + TripodHeight
%   cfg.Mount.Regions, if non-empty, is used as-is (custom layouts).

    if isfield(cfg.Mount, 'Regions') && ~isempty(cfg.Mount.Regions)
        regions = cfg.Mount.Regions;
        return;
    end

    lo = cfg.CamLowerBounds(1:3);
    hi = cfg.CamUpperBounds(1:3);
    footLo = [cfg.Volume(1,1), cfg.Volume(2,1)];
    footHi = [cfg.Volume(1,2), cfg.Volume(2,2)];

    walls = [ ...
        plane('Wall -x',  lo, [lo(1) hi(2) hi(3)]), ...
        plane('Wall +x', [hi(1) lo(2) lo(3)], hi), ...
        plane('Wall -y',  lo, [hi(1) lo(2) hi(3)]), ...
        plane('Wall +y', [lo(1) hi(2) lo(3)], hi), ...
        plane('Ceiling', [lo(1) lo(2) hi(3)], hi)];

    tz = cfg.Mount.TripodHeight;
    r = footRadius(cfg);
    tripod = struct('Name', 'Tripod ring', 'Type', 'ring', ...
        'Lo', [lo(1)+r lo(2)+r tz(1)], 'Hi', [hi(1)-r hi(2)-r tz(2)], ...
        'InnerLo', footLo, 'InnerHi', footHi);

    switch lower(cfg.Mount.Model)
        case 'walls_ceiling_tripod'
            regions = [walls, tripod];
        case 'walls_ceiling'
            regions = walls;
        case 'tripod'
            regions = tripod;
        case 'tripod_cabinets'
            regions = tripodCabinets(cfg, lo, hi, footLo, footHi, r);
        case 'box'
            regions = struct('Name', {}, 'Type', {}, 'Lo', {}, 'Hi', {}, ...
                             'InnerLo', {}, 'InnerHi', {});
        otherwise
            error('mountRegions:UnknownModel', 'Unknown mount model "%s".', cfg.Mount.Model);
    end
end

function regions = tripodCabinets(cfg, lo, hi, footLo, footHi, r)
% Floor band on the x sides (between the footprint and the cabinets), the
% y-end bands (floor or table, wider height range) and the cabinet tops.
% r: tripod foot radius, kept clear of the walls (y ends), the cabinet
% fronts (floor band) and the cabinet-top edges.
    tz = cfg.Mount.TripodHeight;  ez = cfg.Mount.EndHeight;  c = cfg.Mount.Cabinet;
    bx = c.XDepth(1);                         % floor band ends where the cabinets start
    cx = [c.XDepth(1) + r, c.XDepth(2) - r];  cy = [c.YSpan(1) + r, c.YSpan(2) - r];
    if cx(1) > cx(2) || cy(1) > cy(2)
        error('mountRegions:TripodTooWide', ...
            'Tripod foot radius %.2f m does not fit on the %.2f m deep cabinets.', r, diff(c.XDepth));
    end
    regions = [ ...
        box('Floor band -x', [-bx + r  lo(2) + r tz(1)], [footLo(1) hi(2) - r tz(2)]), ...
        box('Floor band +x', [footHi(1) lo(2) + r tz(1)], [bx - r   hi(2) - r tz(2)]), ...
        box('End -y (floor/table)', [-bx lo(2) + r ez(1)], [bx footLo(2) ez(2)]), ...
        box('End +y (floor/table)', [-bx footHi(2) ez(1)], [bx hi(2) - r ez(2)]), ...
        box('Cabinet -x', [-cx(2) cy(1) c.Height + tz(1)], [-cx(1) cy(2) c.Height + tz(2)]), ...
        box('Cabinet +x', [ cx(1) cy(1) c.Height + tz(1)], [ cx(2) cy(2) c.Height + tz(2)])];
end

function r = footRadius(cfg)
    r = 0;
    if isfield(cfg.Mount, 'TripodFootRadius') && ~isempty(cfg.Mount.TripodFootRadius)
        r = cfg.Mount.TripodFootRadius;
    end
end

function r = box(name, a, b)
    r = struct('Name', name, 'Type', 'box', 'Lo', min(a, b), 'Hi', max(a, b), ...
               'InnerLo', [], 'InnerHi', []);
end

function r = plane(name, a, b)
    r = struct('Name', name, 'Type', 'plane', 'Lo', min(a, b), 'Hi', max(a, b), ...
               'InnerLo', [], 'InnerHi', []);
end
