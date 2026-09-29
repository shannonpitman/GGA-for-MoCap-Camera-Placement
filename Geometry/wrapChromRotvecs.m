function chrom = wrapChromRotvecs(chrom)
%WRAPCHROMROTVECS  Apply rotvecWrap to every camera's orientation genes (4-6).
    for s = 0:6:numel(chrom) - 6
        chrom(s+4:s+6) = rotvecWrap(chrom(s+4:s+6));
    end
end
