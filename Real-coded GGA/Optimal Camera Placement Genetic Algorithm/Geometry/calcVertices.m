function vertices = calcVertices(numVisible, visibleIdx, adj, planes) %#ok<INUSD>
% CALCVERTICES  Vertices of the intersection of all visible camera frustums.
%
% A vertex of a convex polytope defined by half-spaces is the intersection
% of three of its bounding planes that also satisfies every other plane.
% This implementation enumerates every triple of distinct planes across
% all visible cameras and keeps the well-conditioned intersections that
% lie inside every frustum.
%
% The hot loops are fully vectorised:
%   * 3x3 plane-triple linear solves are done with Cramer's rule on T
%     triples at once (no per-triple A\b call).
%   * The "is the candidate inside every camera's frustum?" test is a
%     single (P x T) = Nall' * X matmul plus an `any` reduction, where P
%     is the total plane count across visible cameras.
%
% INPUTS
%   numVisible  - number of cameras seeing the target point
%   visibleIdx  - indices into the original `planes` cell of those cameras
%   adj         - kept for API compatibility (not used here)
%   planes      - cell array of (4 x K) plane stacks per camera; columns
%                 are [n; d] for each pyramid surface
%
% OUTPUT
%   vertices    - V x 3 array of intersection points lying inside the
%                 intersection polytope (may be empty if numVisible < 2
%                 or the polytope is degenerate)

    tol     = -1e-9;     % "inside" tolerance for boundary admission
    detThr  = 1e-9;      % "non-parallel" threshold for the 3-plane triple

    if numVisible < 2
        vertices = zeros(0, 3);
        return;
    end

    % --- Concatenate all visible cameras' planes into Nall (3xP), Dall (1xP) ---
    nPerCam = zeros(numVisible, 1);
    for i = 1:numVisible
        nPerCam(i) = size(planes{visibleIdx(i)}, 2);
    end
    P    = sum(nPerCam);
    Nall = zeros(3, P);
    Dall = zeros(1, P);
    col  = 0;
    for i = 1:numVisible
        c = nPerCam(i);
        Nall(:, col+1:col+c) = planes{visibleIdx(i)}(1:3, :);
        Dall(   col+1:col+c) = planes{visibleIdx(i)}(4,   :);
        col = col + c;
    end

    if P < 3
        vertices = zeros(0, 3);
        return;
    end

    % --- Enumerate all triples of distinct planes ---
    triples  = nchoosek(1:P, 3);
    nTriples = size(triples, 1);

    % Pull the three normals and offsets per triple as 3xT / 1xT vectors.
    aIdx = triples(:,1).'; bIdx = triples(:,2).'; cIdx = triples(:,3).';
    n_a = Nall(:, aIdx);  n_b = Nall(:, bIdx);  n_c = Nall(:, cIdx);
    d_a = Dall(aIdx);     d_b = Dall(bIdx);     d_c = Dall(cIdx);

    % --- Vectorised determinant filter ---
    nax = n_a(1,:); nay = n_a(2,:); naz = n_a(3,:);
    nbx = n_b(1,:); nby = n_b(2,:); nbz = n_b(3,:);
    ncx = n_c(1,:); ncy = n_c(2,:); ncz = n_c(3,:);

    % Helpful 2x2 cofactors reused below.
    bcyz = nby.*ncz - nbz.*ncy;
    bcxz = nbx.*ncz - nbz.*ncx;
    bcxy = nbx.*ncy - nby.*ncx;

    detA = nax.*bcyz - nay.*bcxz + naz.*bcxy;

    keep = abs(detA) > detThr;
    if ~any(keep)
        vertices = zeros(0, 3);
        return;
    end

    % Restrict everything to the well-conditioned triples.
    nax = nax(keep); nay = nay(keep); naz = naz(keep);
    nbx = nbx(keep); nby = nby(keep); nbz = nbz(keep);
    ncx = ncx(keep); ncy = ncy(keep); ncz = ncz(keep);
    d_a = d_a(keep); d_b = d_b(keep); d_c = d_c(keep);
    detA = detA(keep);
    bcyz = bcyz(keep); bcxz = bcxz(keep); bcxy = bcxy(keep);

    % --- Vectorised Cramer's rule for x = A \ b across all kept triples ---
    % A_t (3x3) has rows [n_a; n_b; n_c], b_t = [d_a; d_b; d_c].
    %   x = det(A with col 1 replaced by b) / det(A)
    %   y = det(A with col 2 replaced by b) / det(A)
    %   z = det(A with col 3 replaced by b) / det(A)
    detX = d_a.*bcyz                                    ...
         - nay.*(d_b.*ncz - nbz.*d_c)                   ...
         + naz.*(d_b.*ncy - nby.*d_c);

    detY = nax.*(d_b.*ncz - nbz.*d_c)                   ...
         - d_a.*bcxz                                    ...
         + naz.*(nbx.*d_c - d_b.*ncx);

    detZ = nax.*(nby.*d_c - d_b.*ncy)                   ...
         - nay.*(nbx.*d_c - d_b.*ncx)                   ...
         + d_a.*bcxy;

    X = [detX; detY; detZ] ./ detA;     % 3 x T_kept

    % --- Vectorised inside-every-frustum test ---
    % violations(p, t) = n_p . x_t - d_p ; point t is outside if any p violates.
    violations = Nall.' * X - Dall.';   % P x T_kept (broadcast on Dall')
    inside     = ~any(violations < tol, 1);   % 1 x T_kept

    vertices = X(:, inside).';          % V x 3
end
