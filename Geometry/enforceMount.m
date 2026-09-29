function chrom = enforceMount(chrom, regions)
%ENFORCEMOUNT  Project every camera in a chromosome onto a mountable region.
%   Genes 1-3 of each 6-gene camera block are the position; orientation
%   genes are unchanged. No-op when regions is empty.
    if isempty(regions)
        return;
    end
    for s = 0:6:numel(chrom) - 6
        chrom(s+1:s+3) = projectToMount(chrom(s+1:s+3), regions);
    end
end
