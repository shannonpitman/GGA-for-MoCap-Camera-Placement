function d = pairDistances(P)
%PAIRDISTANCES  Euclidean distance between every pair of rows of P.
%   d = pairDistances(P) is a row vector in the same order as pdist(P)
%   (pairs (1,2), (1,3), ..., (n-1,n)). Base MATLAB only, so the GA does
%   not need a Statistics Toolbox licence (the sim machine's licence
%   checkout failed once, 2026-11-26).
    n = size(P, 1);
    [j, i] = find(tril(true(n), -1));   % column-major: i < j, i varies slowest
    d = vecnorm(P(i, :) - P(j, :), 2, 2).';
end
