function chrom = eulerChromToRotvec(chrom)
%EULERCHROMTOROTVEC  Convert a legacy chromosome (XYZ Euler genes, as saved
%   before the exponential-map change) to rotation-vector genes. Positions
%   are unchanged; the rotation of every camera is identical.
    for s = 0:6:numel(chrom) - 6
        chrom(s+4:s+6) = rotmToGenes(eul2rotm(chrom(s+4:s+6), "XYZ"));
    end
end
