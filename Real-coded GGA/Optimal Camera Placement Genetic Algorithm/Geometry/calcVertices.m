function vertices = calcVertices(numVisible, visibleIdx, adj, planes) %#ok<INUSD>
% CALCVERTICES  Vertices of the intersection of all visible camera frustums.
%
% A vertex of a convex polytope defined by half-spaces is the intersection
% of three of its bounding planes that also satisfies every other plane.
% Earlier versions of this function only enumerated triples of the form
% (2 planes from camera I + 1 from camera K) for each unordered camera
% pair, which silently dropped:
%   - the symmetric (1 from I + 2 from K) triples, and
%   - all triples involving three different cameras.
%
% After fixing isInsidePlanes (it now correctly rejects candidates outside
% the polytope), the missing triples started causing degenerate vertex
% sets and NaN uncertainties for some chromosomes. This implementation
% enumerates every triple of distinct planes across all visible cameras
% and keeps the well-conditioned intersections that lie inside every
% camera's frustum.
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

    tol     = -1e-9;
    detThr  = 1e-9;     % positive threshold for "non-parallel" planes

    if numVisible < 2
        vertices = zeros(0, 3);
        return;
    end

    % --- Build a single concatenated plane list across visible cameras ---
    visiblePlanes = cell(numVisible, 1);
    nPerCam       = zeros(numVisible, 1);
    for i = 1:numVisible
        visiblePlanes{i} = planes{visibleIdx(i)};
        nPerCam(i)       = size(visiblePlanes{i}, 2);
    end
    P    = sum(nPerCam);                       % total plane count
    Nall = zeros(3, P);
    Dall = zeros(1, P);
    col  = 0;
    for i = 1:numVisible
        c = nPerCam(i);
        Nall(:, col+1:col+c) = visiblePlanes{i}(1:3, :);
        Dall(   col+1:col+c) = visiblePlanes{i}(4,   :);
        col = col + c;
    end

    if P < 3
        vertices = zeros(0, 3);
        return;
    end

    % --- Enumerate all triples of distinct planes ---
    triples  = nchoosek(1:P, 3);
    nTriples = size(triples, 1);

    % --- Vectorised determinant filter (skip parallel/ill-conditioned) ---
    % A_t = [n_a; n_b; n_c]^T (3x3); det via the cofactor expansion below.
    aIdx = triples(:,1); bIdx = triples(:,2); cIdx = triples(:,3);
    n_a = Nall(:, aIdx);  % 3 x T
    n_b = Nall(:, bIdx);
    n_c = Nall(:, cIdx);
    % det([n_a n_b n_c]) for each column-triple
    detTriples = n_a(1,:) .* (n_b(2,:).*n_c(3,:) - n_b(3,:).*n_c(2,:)) ...
               - n_a(2,:) .* (n_b(1,:).*n_c(3,:) - n_b(3,:).*n_c(1,:)) ...
               + n_a(3,:) .* (n_b(1,:).*n_c(2,:) - n_b(2,:).*n_c(1,:));
    keep = abs(detTriples) > detThr;
    if ~any(keep)
        vertices = zeros(0, 3);
        return;
    end

    triplesK = triples(keep, :);
    out   = zeros(size(triplesK,1), 3);
    nKept = 0;
    for t = 1:size(triplesK, 1)
        idxT = triplesK(t, :);
        A = Nall(:, idxT).';     % 3x3 plane normals as rows
        b = Dall(idxT).';        % 3x1
        x = A \ b;               % candidate vertex

        if isInsidePlanes(x, visiblePlanes, tol, numVisible)
            nKept           = nKept + 1;
            out(nKept, :)   = x.';
        end
    end

    vertices = out(1:nKept, :);
end
