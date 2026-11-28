function [chrom, nViol] = separateTripods(chrom, regions, minSpacing, maxIter)
%SEPARATETRIPODS  Push cameras apart until every pair is minSpacing apart in x-y.
%
%   [chrom, nViol] = separateTripods(chrom, regions, minSpacing)
%
%   Every camera stands on its own tripod, so two tripod centres closer than
%   minSpacing [m] (horizontal distance) cannot both be built. Each
%   iteration moves both cameras of every too-close pair apart along the
%   line joining them (half the shortfall each) and projects them back onto
%   the mount regions (projectToMount). Orientation genes are unchanged.
%   nViol is the number of pairs still too close after maxIter iterations
%   (default 100); 0 when the repair succeeded. No-op for minSpacing <= 0.

    if nargin < 4 || isempty(maxIter), maxIter = 100; end
    nViol = 0;
    if isempty(minSpacing) || minSpacing <= 0
        return;
    end
    n = numel(chrom) / 6;
    P = reshape(chrom, 6, n).';
    tol = 1e-6;
    for it = 1:maxIter
        move = zeros(n, 2);
        nViol = 0;
        for i = 1:n-1
            for j = i+1:n
                v = P(i, 1:2) - P(j, 1:2);
                d = norm(v);
                if d >= minSpacing - tol, continue; end
                nViol = nViol + 1;
                if d < 1e-9
                    a = 2*pi*(i + 0.618*j);        % coincident: fixed, pair-specific direction
                    u = [cos(a) sin(a)];
                else
                    u = v / d;
                end
                push = (minSpacing - d)/2 + 1e-3;
                move(i, :) = move(i, :) + push*u;
                move(j, :) = move(j, :) - push*u;
            end
        end
        if nViol == 0, break; end
        P(:, 1:2) = P(:, 1:2) + move;
        for k = 1:n
            P(k, 1:3) = projectToMount(P(k, 1:3), regions);
        end
    end
    if nViol > 0                                    % count after the last move
        nViol = nnz(pairDistances(P(:, 1:2)) < minSpacing - tol);
    end
    chrom = reshape(P.', 1, []);
end
