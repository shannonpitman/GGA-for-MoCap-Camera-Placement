function chrom = enforceMount(chrom, regions, minSpacing, minElevDeg)
%ENFORCEMOUNT  Make every camera in a chromosome physically mountable.
%   chrom = enforceMount(chrom, regions)
%   chrom = enforceMount(chrom, regions, minSpacing)
%   chrom = enforceMount(chrom, regions, minSpacing, minElevDeg)
%
%   Genes 1-3 of each 6-gene camera block (position) are projected onto
%   the mountable regions (no-op when regions is empty). With
%   minSpacing > 0 (cfg.Mount.MinTripodSpacing), cameras closer than that
%   in x-y are then pushed apart (separateTripods). With minElevDeg
%   (cfg.Mount.MinElevation, e.g. 0 = level or up), cameras pointing lower
%   are tilted up to the limit (limitElevation); otherwise the orientation
%   genes are unchanged.
    if nargin >= 4
        chrom = limitElevation(chrom, minElevDeg);
    end
    if isempty(regions)
        return;
    end
    for s = 0:6:numel(chrom) - 6
        chrom(s+1:s+3) = projectToMount(chrom(s+1:s+3), regions);
    end
    if nargin >= 3 && ~isempty(minSpacing) && minSpacing > 0
        chrom = separateTripods(chrom, regions, minSpacing);
    end
end
