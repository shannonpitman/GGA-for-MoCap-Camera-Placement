function eulDeg = chromToEulerDeg(chrom)
%CHROMTOEULERDEG  numCams x 3 XYZ Euler angles [deg] for display/reporting.
%   R = Rx(a) Ry(b) Rz(g), the convention used by the pre-2026-09 results.
    n = numel(chrom) / 6;
    eulDeg = zeros(n, 3);
    for c = 1:n
        eulDeg(c, :) = rad2deg(rotm2eul(genesToRotm(chrom((c-1)*6+4:c*6)), "XYZ"));
    end
end
