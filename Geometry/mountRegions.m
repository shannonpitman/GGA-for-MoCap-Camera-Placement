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
%     'box'                   no restriction beyond the search box (legacy)
%   cfg.Mount.TripodHeight = [zMin zMax] sets the tripod head height range.
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
    tripod = struct('Name', 'Tripod ring', 'Type', 'ring', ...
        'Lo', [lo(1) lo(2) tz(1)], 'Hi', [hi(1) hi(2) tz(2)], ...
        'InnerLo', footLo, 'InnerHi', footHi);

    switch lower(cfg.Mount.Model)
        case 'walls_ceiling_tripod'
            regions = [walls, tripod];
        case 'walls_ceiling'
            regions = walls;
        case 'tripod'
            regions = tripod;
        case 'box'
            regions = struct('Name', {}, 'Type', {}, 'Lo', {}, 'Hi', {}, ...
                             'InnerLo', {}, 'InnerHi', {});
        otherwise
            error('mountRegions:UnknownModel', 'Unknown mount model "%s".', cfg.Mount.Model);
    end
end

function r = plane(name, a, b)
    r = struct('Name', name, 'Type', 'plane', 'Lo', min(a, b), 'Hi', max(a, b), ...
               'InnerLo', [], 'InnerHi', []);
end
