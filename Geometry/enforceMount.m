function chrom = enforceMount(chrom, regions, minSpacing)
%ENFORCEMOUNT  Project every camera in a chromosome onto a mountable region.
%   chrom = enforceMount(chrom, regions)
%   chrom = enforceMount(chrom, regions, minSpacing)
%
%   Genes 1-3 of each 6-gene camera block are the position; orientation
%   genes are unchanged. No-op when regions is empty. With minSpacing > 0
%   (cfg.Mount.MinTripodSpacing), cameras closer than that in x-y are then
%   pushed apart (separateTripods).
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
