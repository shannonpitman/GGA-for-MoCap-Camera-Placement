function feas = feasibilityMap(opts)
%FEASIBILITYMAP  Config-independent feasibility of >=2 triangulable coverage.
%
%   feas = feasibilityMap()        % defaults: 6 cams, full 3D flight volume
%   feas = feasibilityMap(opts)    % override defaults via a struct
%
% PURPOSE
%   Answers the question you must settle BEFORE wiring a >=2 coverage
%   requirement into the GA: "Which target points COULD ever be covered by
%   >=2 triangulable cameras, given only the camera-mounting region --
%   independent of any particular configuration?"
%
%   A point that fails this test can NEVER be double-covered by ANY 6-camera
%   arrangement, so a HARD >=2 constraint would make every config infeasible
%   and break the GA. If any such dead-zone points exist you must instead use
%   a soft penalty (or exclude those points from the constraint).
%
% METHOD  (deliberately a best-case UPPER BOUND)
%   * Candidate camera positions are sampled over the mounting surfaces
%     (the 4 walls + ceiling of the position-bound box). Surface positions
%     give the widest possible triangulation baselines, so this is the most
%     optimistic admissible candidate set. Interior positions would only
%     shrink baselines (worse), so restricting to surfaces is correct for an
%     upper bound.
%   * Each candidate is assumed AIMED at the point (best case), so the point
%     sits at the principal point and field-of-view is always satisfied --
%     only effective RANGE and the parallax ANGLE bind. This makes the result
%     a NECESSARY condition: real 6-camera configs (fixed orientations, the
%     same 6 cameras shared across every point) will always do worse, never
%     better than this map.
%   * A point is "double-coverable" if there exist >=2 candidate positions
%     within range whose camera->point directions subtend a convergence
%     angle inside [minTriangAngle, maxTriangAngle] -- the SAME angle test as
%     checkTriangulability.m (angle between the two camera->point unit
%     vectors).
%
% NOTE ON THE COVERAGE DEFINITION
%   This map uses >=2 TRIANGULABLE views (the meaningful mocap bar you
%   selected). The current CF2 hard penalty (calculatePointOcclusion.m) fires
%   on <2 VISIBLE views, which is strictly weaker. The gap between the
%   "feasibleVisible" and "feasible" (triangulable) fields below tells you how
%   much stricter the triangulable requirement is for this volume.
%
% RANGE REGIMES (odd cams = narrow 16 m, even cams = wide 9 m)
%   * feasible       : narrow-lens range (16 m) allowed at every candidate.
%                      Failures here are HARD geometric dead zones.
%   * feasibleWide   : wide-lens range (9 m). Points feasible only in the
%                      narrow pass are "narrowDependent" -- coverable only if
%                      >=2 of the 3 long-range cameras are spent on them.
%
% OUTPUT struct feas has one row per target point (feas.points, N x 3):
%   feas.feasible        N x 1 logical  HARD >=2 triangulable feasibility (16 m)
%   feas.feasibleWide    N x 1 logical  >=2 triangulable within wide range (9 m)
%   feas.feasibleVisible N x 1 logical  >=2 merely-visible (what CF2 checks)
%   feas.narrowDependent N x 1 logical  feasible only via long-range cameras
%   feas.nInRangeNarrow  N x 1  # candidate mounts within 16 m
%   feas.nInRangeWide    N x 1  # candidate mounts within 9 m
%   feas.maxParallax     N x 1  largest parallax angle available (16 m) [deg]
%
% Requires the project functions on path (setupHardwareSpecs, setupCostParams,
% generateTargetSpace); addProjectPaths is called if available.

    if nargin < 1, opts = struct; end
    if exist('addProjectPaths','file'), addProjectPaths(); end

    % ---- defaults (match runCameraOptimiser.m) --------------------------
    def = struct( ...
        'numCams',      6, ...
        'targetType',   1, ...             % 1 = full UAV flight volume (full 3D)
        'targetMode',   1, ...             % 1 = uniform grid
        'spacing',      1, ...             % target grid spacing [m]
        'volume',       [-4 4; -4 4; 0 4], ...   % MS.G flight envelope [m]
        'camLower',     [-5 -4.5 0.0], ...       % camera position lower bound
        'camUpper',     [ 5  4.5 4.8], ...       % camera position upper bound
        'candStep',     0.5, ...           % mounting-surface sample step [m]
        'includeFloor', false, ...         % allow floor-mounted candidates?
        'ugvSlabZ',     0.5, ...           % z<=this is the UGV floor slab [m]
        'saveFig',      true, ...
        'saveMat',      true, ...
        'outDir',       fullfile(pwd,'Results','Feasibility'));
    fn = fieldnames(def);
    for k = 1:numel(fn)
        if ~isfield(opts,fn{k}) || isempty(opts.(fn{k}))
            opts.(fn{k}) = def.(fn{k});
        end
    end

    % ---- constants from the project's OWN setup (single source of truth)-
    specs = setupHardwareSpecs(opts.numCams);
    specs.TargetType   = opts.targetType;
    specs.TargetMode   = opts.targetMode;
    specs.Target       = generateTargetSpace(opts.volume, opts.targetMode, opts.spacing);
    specs.NumPoints    = size(specs.Target,1);
    specs.spacing      = opts.spacing;
    specs.UseNormTable = false;            % feasibility needs no normTable
    specs = setupCostParams(specs);

    minAng  = specs.PreComputed.minTriangAngle;      % 40 deg
    maxAng  = specs.PreComputed.maxTriangAngle;      % 140 deg
    rNarrow = specs.PreComputed.maxCameraRange;      % 16 m
    rWide   = specs.PreComputed.maxCameraRangeWide;  % 9 m

    P = specs.Target;                      % N x 3 target points
    N = size(P,1);

    % ---- candidate camera positions on the mounting surfaces ------------
    C = mountingCandidates(opts.camLower, opts.camUpper, opts.candStep, opts.includeFloor);
    M = size(C,1);

    % ---- per-point feasibility ------------------------------------------
    feasible       = false(N,1);   % >=2 triangulable, narrow range (HARD bound)
    feasibleWide   = false(N,1);   % >=2 triangulable, wide range
    feasibleVis    = false(N,1);   % >=2 visible (what CF2 currently penalises)
    nInRangeNarrow = zeros(N,1);
    nInRangeWide   = zeros(N,1);
    maxParallax    = zeros(N,1);

    for p = 1:N
        v   = P(p,:) - C;                    % M x 3 : camera -> point vectors
        rng = sqrt(sum(v.^2,2));             % M x 1 : range to point
        u   = v ./ max(rng, eps);            % unit camera->point directions

        inN = rng <= rNarrow & rng > 0;      % within narrow-lens range
        inW = rng <= rWide   & rng > 0;      % within wide-lens range
        nInRangeNarrow(p) = sum(inN);
        nInRangeWide(p)   = sum(inW);

        feasibleVis(p) = nInRangeNarrow(p) >= 2;   % >=2 can see it at all

        [feasible(p),  maxParallax(p)] = hasTriangPair(u(inN,:), minAng, maxAng);
        feasibleWide(p)                = hasTriangPair(u(inW,:), minAng, maxAng);
    end

    narrowDependent = feasible & ~feasibleWide;

    % ---- pack results ---------------------------------------------------
    feas.points          = P;
    feas.feasible        = feasible;
    feas.feasibleWide    = feasibleWide;
    feas.feasibleVisible = feasibleVis;
    feas.narrowDependent = narrowDependent;
    feas.nInRangeNarrow  = nInRangeNarrow;
    feas.nInRangeWide    = nInRangeWide;
    feas.maxParallax     = maxParallax;
    feas.candidates      = C;
    feas.params = struct('minAng',minAng,'maxAng',maxAng, ...
        'rNarrow',rNarrow,'rWide',rWide,'candStep',opts.candStep, ...
        'spacing',opts.spacing,'numCams',opts.numCams);

    % ---- summary --------------------------------------------------------
    nInf = sum(~feasible);
    fprintf('\n=== >=2 Triangulable Feasibility (%d cams, best-case upper bound) ===\n', opts.numCams);
    fprintf('Target points               : %d\n', N);
    fprintf('Candidate mount positions   : %d  (surface step %.2f m)\n', M, opts.candStep);
    fprintf('Triangulation angle window  : [%g, %g] deg\n', minAng, maxAng);
    fprintf('Ranges (narrow / wide)      : %g / %g m\n', rNarrow, rWide);
    fprintf('---------------------------------------------------------------\n');
    fprintf('>=2 visible (CF2 bar)       : %d (%.1f%%) feasible\n', ...
        sum(feasibleVis), 100*sum(feasibleVis)/N);
    fprintf('>=2 triangulable (16 m)     : %d (%.1f%%) feasible\n', ...
        sum(feasible), 100*sum(feasible)/N);
    fprintf('HARD-infeasible points      : %d (%.1f%%)  <- NO config can double-cover these\n', ...
        nInf, 100*nInf/N);
    fprintf('Narrow-lens-dependent points: %d (%.1f%%)  <- need >=2 of the 3 long-range cams\n', ...
        sum(narrowDependent), 100*sum(narrowDependent)/N);
    if nInf == 0
        fprintf('=> A HARD >=2 triangulable constraint is ADMISSIBLE for this volume.\n');
    else
        fprintf('=> HARD constraint NOT admissible; use a soft penalty or exclude the %d dead-zone points.\n', nInf);
    end

    slab = P(:,3) <= opts.ugvSlabZ;
    if any(slab)
        fprintf('UGV floor slab (z<=%.2fm)   : %d pts, %d hard-infeasible (%.1f%%)\n', ...
            opts.ugvSlabZ, sum(slab), sum(~feasible & slab), ...
            100*sum(~feasible & slab)/max(1,sum(slab)));
    end
    fprintf('===============================================================\n');

    % ---- figure + save --------------------------------------------------
    if opts.saveFig,  plotFeasibility(feas, opts); end
    if opts.saveMat
        if ~exist(opts.outDir,'dir'), mkdir(opts.outDir); end
        matName = fullfile(opts.outDir, ...
            sprintf('Feasibility_%dCam_sp%03.0fcm.mat', opts.numCams, 100*opts.spacing));
        save(matName, 'feas', 'opts');
        fprintf('Data saved  : %s\n', matName);
    end
end

% ======================================================================= %
function [ok, maxAngle] = hasTriangPair(U, minAng, maxAng)
%HASTRIANGPAIR  True if any pair of unit dirs U (K x 3) has convergence angle
%   in [minAng, maxAng]. Convergence angle = acosd(dot(u_i,u_j)), matching
%   checkTriangulability.m. maxAngle returns the largest available pair angle.
    ok = false; maxAngle = 0;
    K = size(U,1);
    if K < 2, return; end
    G = U * U.';                       % K x K dot products of unit vectors
    G = min(1, max(-1, G));            % clamp for acosd
    A = acosd(G);                      % pairwise convergence angles [deg]
    A(1:K+1:end) = NaN;                % drop self-pairs on the diagonal
    maxAngle = max(A(:));
    ok = any(A(:) >= minAng & A(:) <= maxAng, 'all');
end

% ======================================================================= %
function C = mountingCandidates(lo, hi, step, includeFloor)
%MOUNTINGCANDIDATES  Sample the 4 walls + ceiling (optionally floor) of the
%   camera-position box [lo; hi] at the given step. Returns unique M x 3 pts.
    xs = axisSamples(lo(1), hi(1), step);
    ys = axisSamples(lo(2), hi(2), step);
    zs = axisSamples(lo(3), hi(3), step);

    C = zeros(0,3);
    % Walls X = lo(1) and hi(1)  (vary Y,Z)
    [Yg,Zg] = meshgrid(ys, zs);
    C = [C; lo(1)*ones(numel(Yg),1), Yg(:), Zg(:)];
    C = [C; hi(1)*ones(numel(Yg),1), Yg(:), Zg(:)];
    % Walls Y = lo(2) and hi(2)  (vary X,Z)
    [Xg,Zg] = meshgrid(xs, zs);
    C = [C; Xg(:), lo(2)*ones(numel(Xg),1), Zg(:)];
    C = [C; Xg(:), hi(2)*ones(numel(Xg),1), Zg(:)];
    % Ceiling Z = hi(3)  (vary X,Y)
    [Xg,Yg] = meshgrid(xs, ys);
    C = [C; Xg(:), Yg(:), hi(3)*ones(numel(Xg),1)];
    if includeFloor
        C = [C; Xg(:), Yg(:), lo(3)*ones(numel(Xg),1)];
    end
    C = unique(C, 'rows');
end

% ======================================================================= %
function s = axisSamples(a, b, step)
    s = a:step:b;
    if isempty(s) || s(end) ~= b, s = [s, b]; end
end

% ======================================================================= %
function plotFeasibility(feas, opts)
    P  = feas.points;
    bad = ~feas.feasible;
    nd  = feas.narrowDependent;
    ok  = feas.feasible & ~nd;

    fig = figure('Color','w','Position',[100 100 960 720]);
    hold on; grid on; axis equal; box on;
    scatter3(feas.candidates(:,1), feas.candidates(:,2), feas.candidates(:,3), ...
        6, [0.6 0.6 0.85], '.');
    if any(ok)
        scatter3(P(ok,1),  P(ok,2),  P(ok,3),  28, [0.20 0.70 0.25], 'filled', ...
            'MarkerFaceAlpha', 0.45);
    end
    if any(nd)
        scatter3(P(nd,1),  P(nd,2),  P(nd,3),  46, [0.95 0.60 0.10], 'filled');
    end
    if any(bad)
        scatter3(P(bad,1), P(bad,2), P(bad,3), 70, [0.85 0.10 0.10], 'filled');
    end
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    title(sprintf('>=2 triangulable feasibility  |  %d cams  |  %d hard dead-zone pts (%.1f%%)', ...
        opts.numCams, sum(bad), 100*sum(bad)/numel(bad)));
    legend({'mount candidates','feasible','narrow-lens dependent','HARD infeasible'}, ...
        'Location','bestoutside');
    view(45,25);

    if opts.saveFig
        if ~exist(opts.outDir,'dir'), mkdir(opts.outDir); end
        figName = fullfile(opts.outDir, ...
            sprintf('FeasibilityMap_%dCam_sp%03.0fcm.pdf', opts.numCams, 100*opts.spacing));
        try
            exportgraphics(fig, figName, 'ContentType','vector');
            fprintf('Figure saved: %s\n', figName);
        catch
            saveas(fig, strrep(figName,'.pdf','.fig'));
            fprintf('Figure saved (.fig fallback): %s\n', strrep(figName,'.pdf','.fig'));
        end
    end
end
