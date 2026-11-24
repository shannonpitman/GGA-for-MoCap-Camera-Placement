function [chrom, info] = legaliseRig(chrom0, specs, preset)
%LEGALISERIG  Move an existing rig onto a preset's tripod-footprint mounts and re-cost it.
%
%   [chrom, info] = legaliseRig(chrom0, specs, preset)
%
%   Post-processing for rigs optimised without the tripod footprint (e.g.
%   openmv_cab_2m): each camera is projected onto the mount regions of
%   preset (e.g. 'openmv_cab_2m_fp': legs clear of walls and cabinet edges)
%   and tripods closer than Mount.MinTripodSpacing are pushed apart
%   (enforceMount / separateTripods). Cameras pointing below
%   Mount.MinElevation are tilted up to it (limitElevation); other
%   orientations are kept as they are.
%   specs is the run's saveData.Specifications, used to re-cost the rig on
%   the same target grid.
%
%   info fields:
%     Moved          n x 1 horizontal+vertical move of each camera [m]
%     MaxMove        largest move [m]
%     Violations     tripod pairs still too close (0 = repaired)
%     RawUnc_mm      [before after] resolution uncertainty
%     RawOcc_deg     [before after] dynamic occlusion
%     TwoPlusPct     [before after] % of target points seen by 2+ cameras

    cfg = runConfig(preset);
    regions = mountRegions(cfg);
    chrom = enforceMount(chrom0, regions, cfg.Mount.MinTripodSpacing, cfg.Mount.MinElevation);
    P0 = reshape(chrom0, 6, []).';  P1 = reshape(chrom, 6, []).';
    nViol = nnz(pdist(P1(:, 1:2)) < cfg.Mount.MinTripodSpacing - 1e-6);
    info.Moved = vecnorm(P1(:, 1:3) - P0(:, 1:3), 2, 2);
    info.MaxMove = max(info.Moved);
    info.Violations = nViol;

    b0 = costBreakdown(chrom0, specs);  b1 = costBreakdown(chrom, specs);
    [~, s0] = perTargetCoverage(chrom0, specs);  [~, s1] = perTargetCoverage(chrom, specs);
    info.RawUnc_mm  = 1000*[b0.RawUnc b1.RawUnc];
    info.RawOcc_deg = [b0.RawOcc b1.RawOcc];
    info.TwoPlusPct = [s0.twoPlusPct s1.twoPlusPct];
end
