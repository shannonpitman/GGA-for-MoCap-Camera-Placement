function specs = setupCostParams(specs)

% Pre-compute data for resolution uncertainty
specs.PreComputed.adjacentSurfaces = [1 2; 2 3; 3 4; 4 1];
specs.PreComputed.du = 0.5; %pixels
specs.PreComputed.dv = 0.5; %pixels
specs.PreComputed.penaltyUncertainty = 1.5; % penalty for unseen points
specs.PreComputed.w2 = 0.2; %weight for ellipsoid optimization

% Pre-compute data for dynamic occlusion
specs.PreComputed.minTriangAngle = 40; % min convergence angle [degrees]
specs.PreComputed.maxTriangAngle = 140; % max convergence angle [degrees]
specs.PreComputed.maxCameraRange = specs.Range; %m narrowFOV effective range
specs.PreComputed.maxCameraRangeWide = specs.RangeWide; %m wideFOV effective range

% Pre-compute target space in homogeneous coordinates for batch processing
specs.TargetHomogeneous = [specs.Target'; ones(1, specs.NumPoints)];

% Normalisation constants + utopia offsets for the combined (CF3) cost.
% Each component is min-max scaled as (raw - utopia)/norm, where
% norm = nadir - utopia, so utopia -> 0 and nadir -> 1.
%
% Defaults below are the fallback used when no calibrated normTable entry
% exists for this instance (utopia = 0 reproduces the old raw/norm scaling).
% If a calibrated entry IS found, all four values are overwritten from
% Results/normTable.mat (see buildNormalisationSchedule.m / getNormConstants.m).
specs.PreComputed.uncertNorm   = 0.08;
specs.PreComputed.occlNorm     = 100;
specs.PreComputed.uncertUtopia = 0;
specs.PreComputed.occlUtopia   = 0;
specs.PreComputed.normSource   = 'default';

% Look up calibrated constants unless explicitly disabled (the schedule
% builder sets specs.UseNormTable = false so it never depends on its own
% output). Requires specs.TargetType/TargetMode/Cams/spacing to be set.
useTable = ~isfield(specs, 'UseNormTable') || specs.UseNormTable;
if useTable
    try
        nc = getNormConstants(specs);
        specs.PreComputed.uncertNorm   = nc.uncertNorm;
        specs.PreComputed.occlNorm     = nc.occlNorm;
        specs.PreComputed.uncertUtopia = nc.utopiaUnc;
        specs.PreComputed.occlUtopia   = nc.utopiaOcc;
        specs.PreComputed.normSource   = 'normTable';
    catch ME
        % No file / no matching instance -> keep defaults. Harmless for
        % CF1/CF2 runs (they ignore these), but a CF3 cost on the defaults is
        % on a different scale from every tabled result, so say so (once per
        % instance per session).
        persistent warned
        if isempty(warned), warned = containers.Map(); end
        key = sprintf('%d_%d_%d_%.3f', specs.TargetType, specs.TargetMode, specs.Cams, specs.spacing);
        if ~isKey(warned, key)
            warning('setupCostParams:NoNormEntry', ...
                ['No normTable entry for TT%d GM%d %dC sp%.2f (%s): CF3 uses the ' ...
                 'default constants, not comparable with tabled results.'], ...
                specs.TargetType, specs.TargetMode, specs.Cams, specs.spacing, ME.message);
            warned(key) = true;
        end
    end
end