function d = layoutDistance(chromA, chromB)
%LAYOUTDISTANCE  Mean camera-position distance between two layouts [m].
%
%   d = layoutDistance(chromA, chromB)
%
%   Camera slot order in a chromosome is arbitrary (two layouts with the same
%   cameras in a different order are the same rig), so cameras are first
%   paired one-to-one to minimise total distance (matchpairs, Hungarian
%   assignment), then the mean paired distance is returned. Orientation is
%   ignored: this is a distance between where the cameras are.
    PA = reshape(chromA, 6, []).';  PA = PA(:, 1:3);
    PB = reshape(chromB, 6, []).';  PB = PB(:, 1:3);
    n = size(PA, 1);
    C = zeros(n);
    for i = 1:n
        C(i, :) = vecnorm(PB - PA(i, :), 2, 2).';
    end
    M = matchpairs(C, 1e6);            % every camera matched (unmatched cost huge)
    d = mean(C(sub2ind([n n], M(:, 1), M(:, 2))));
end
