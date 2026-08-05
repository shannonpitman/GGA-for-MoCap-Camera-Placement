function [J, Junc, Jocc] = cf3Terms(costUnc, costOcc, specs)
%CF3TERMS  Weighted, utopia-shifted CF3 terms from raw component costs.
%
% Single source of truth for how the combined (CF3) cost is assembled from
% the raw resolution-uncertainty and dynamic-occlusion values, so the live
% cost function and every downstream re-computation (reportCostBreakdown,
% evaluateOptiTrackCost, spacingSensitivityCore, saveOptiTrackAsRun, ...)
% stay on exactly the same scale.
%
% Each component is min-max scaled as (raw - utopia)/norm, where
% norm = nadir - utopia, so utopia -> 0 and nadir -> 1, then weighted:
%
%   Junc = WeightUncertainty * (costUnc - uncertUtopia)/uncertNorm
%   Jocc = WeightOcclusion   * (costOcc - occlUtopia)  /occlNorm
%   J    = Junc + Jocc
%
% INPUTS:
%   costUnc, costOcc - raw resUncertainty / dynamicOcclusion values.
%   specs            - must carry WeightUncertainty, WeightOcclusion and
%                      PreComputed.uncertNorm / occlNorm. The utopia offsets
%                      PreComputed.uncertUtopia / occlUtopia are optional;
%                      they default to 0 for legacy specs saved before the
%                      utopia-shift convention (reproducing raw/norm scaling).
%
% OUTPUTS:
%   J    - combined weighted CF3 cost.
%   Junc - weighted resolution term (the J_uncertainty part of J).
%   Jocc - weighted occlusion  term (the J_occlusion  part of J).

    pc = specs.PreComputed;

    if isfield(pc, 'uncertUtopia'), uncertUtopia = pc.uncertUtopia; else, uncertUtopia = 0; end
    if isfield(pc, 'occlUtopia'),   occlUtopia   = pc.occlUtopia;   else, occlUtopia   = 0; end

    Junc = specs.WeightUncertainty * ((costUnc - uncertUtopia) / pc.uncertNorm);
    Jocc = specs.WeightOcclusion   * ((costOcc - occlUtopia)   / pc.occlNorm);
    J    = Junc + Jocc;
end
