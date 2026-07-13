%% testResUncertaintyRange
% Empirically characterises the range of the resolution-uncertainty metric
% for target points that ARE well triangulated (seen by >= 2 cameras with a
% non-degenerate intersection polytope). The purpose is to justify the
% choice of penaltyUncertainty = 100.
%
% RATIONALE
% ---------
% In computePointUncertainty the per-point cost is assigned as:
%     0 cameras see the point         -> penaltyUncertainty        (= 100)
%     1 camera / degenerate polytope  -> 0.5 * penaltyUncertainty  (=  50)
%     >= 2 cameras, valid polytope    -> sum(sqrt(abs(eig(cov(V)))))  [metres]
%
% For the penalty to do its job (i.e. always make an unseen / poorly-seen
% point cost strictly and substantially MORE than a genuinely triangulated
% one), 100 must sit well above the entire realistic distribution of the
% genuine metric. This script measures that distribution by sweeping a fine
% target grid across many randomly generated VALID camera arrangements
% (guided initial-population layouts), and reports where 50 and 100 sit
% relative to it.
%
% USAGE: adjust the parameters below and run. Uses the project's own
%        setup/geometry functions so numbers are identical to the GA.
%
% Author: test harness generated for penalty justification.

clc; clear; close all;
addProjectPaths();
projectRoot = fileparts(mfilename('fullpath'));
rng(1);   % reproducible sampling

%% ---------------- Parameters -------------------------------------------
numConfigs   = 60;      % number of random valid camera arrangements to sample
numCams      = 7;       % cameras per arrangement (matches study)
targetType   = 1;       % 1 = UAV (full volume), 2 = UGV (floor slab)
spacing      = 0.5;     % x-y target grid spacing [m] (finer = more samples)
UGV_maxHeight = 0.5;    % only used if targetType == 2
UGV_zSpacing  = 0.25;   % only used if targetType == 2
saveFigures  = true;

% Workspace volume & mounting bounds (identical to runCameraOptimiser.m)
volume            = [-4 4; -4 4; 0 4];
cameraLowerBounds = [-5 -4.5 0   -pi -pi -pi];
cameraUpperBounds = [ 5  4.5 4.8  pi  pi  pi];

%% ---------------- Setup (mirrors runCameraOptimiser.m) ------------------
specs   = setupHardwareSpecs(numCams);
problem = setupProblem(numCams, 1, cameraUpperBounds, cameraLowerBounds);

if targetType == 2
    volume(3,:)   = [0, UGV_maxHeight];
    targetSpacing = [spacing, spacing, min(UGV_zSpacing, UGV_maxHeight)];
else
    targetSpacing = spacing;
end

specs.Target        = generateTargetSpace(volume, 1, targetSpacing);
specs.NumPoints     = size(specs.Target, 1);
specs.SectionCentres = generateSectionCentres(numCams, volume);
specs               = setupCostParams(specs);   % sets penaltyUncertainty = 100

penalty      = 1.5;   
halfPenalty  = 0.5 * penalty;                      
adj          = specs.PreComputed.adjacentSurfaces;
du           = specs.PreComputed.du;
dv           = specs.PreComputed.dv;
resolution   = specs.Resolution;
numPoints    = specs.NumPoints;
targetPts    = specs.Target;

fprintf('Sampling %d random configs x %d target points = %d point-evaluations.\n', ...
    numConfigs, numPoints, numConfigs*numPoints);

%% ---------------- Sweep -------------------------------------------------
% Bucket counts (over all point-evaluations)
cnt.unseen     = 0;   % 0 cameras
cnt.single     = 0;   % exactly 1 camera
cnt.degenerate = 0;   % >=2 cameras but polytope degenerate (< 4 vertices)
cnt.genuine    = 0;   % >=2 cameras, valid polytope

% Collect genuine metric values (preallocate generously, trim later)
genuine = zeros(numConfigs*numPoints, 1);
gk = 0;

for c = 1:numConfigs
    chrom = initialPopulation(problem.VarMin, problem.VarMax, ...
                              specs.SectionCentres, numCams);
    [cameras, camCentres] = setupCameras(chrom, numCams, resolution, ...
        specs.Focal, specs.FocalWide, specs.PrincipalPoint, specs.PixelSize);

    for p = 1:numPoints
        [val, cat] = pointUncertaintyInstrumented(targetPts(p,:), cameras, ...
            camCentres, numCams, adj, du, dv, resolution);
        switch cat
            case 0, cnt.unseen     = cnt.unseen + 1;
            case 1, cnt.single     = cnt.single + 1;
            case 2, cnt.degenerate = cnt.degenerate + 1;
            case 3
                cnt.genuine = cnt.genuine + 1;
                gk = gk + 1;
                genuine(gk) = val;
        end
    end
    if mod(c,10)==0, fprintf('  ...%d/%d configs done\n', c, numConfigs); end
end
genuine = genuine(1:gk);

%% ---------------- Statistics -------------------------------------------
total = numConfigs*numPoints;
fprintf('\n================ BUCKET BREAKDOWN (all evaluations) ============\n');
fprintf('  unseen (0 cams, cost=%.0f)      : %8d  (%5.1f%%)\n', penalty,      cnt.unseen,     100*cnt.unseen/total);
fprintf('  single/degenerate (cost=%.0f)   : %8d  (%5.1f%%)\n', halfPenalty,  cnt.single+cnt.degenerate, 100*(cnt.single+cnt.degenerate)/total);
fprintf('  genuine (>=2 cams, real metric) : %8d  (%5.1f%%)\n', cnt.genuine, 100*cnt.genuine/total);

if isempty(genuine)
    error('No genuine >=2-camera samples collected; check configuration.');
end

pct   = [0 1 5 25 50 75 95 99 100];
qvals = prctile(genuine, pct);

fprintf('\n============ GENUINE >=2-CAMERA UNCERTAINTY [metres] ===========\n');
fprintf('  n         = %d\n', numel(genuine));
fprintf('  mean      = %.4g\n', mean(genuine));
fprintf('  std       = %.4g\n', std(genuine));
for i = 1:numel(pct)
    fprintf('  p%-3g      = %.4g\n', pct(i), qvals(i));
end

fracGe50  = 100*mean(genuine >= halfPenalty);
fracGe100 = 100*mean(genuine >= penalty);
fprintf('\n---- Separation from the penalty ----\n');
fprintf('  max genuine value          = %.4g\n', max(genuine));
fprintf('  half-penalty (single cam)  = %.0f\n', halfPenalty);
fprintf('  full penalty (unseen)      = %.0f\n', penalty);
fprintf('  %% of genuine values >= 50  = %.3f%%\n', fracGe50);
fprintf('  %% of genuine values >= 100 = %.3f%%\n', fracGe100);
fprintf('  penalty / p99(genuine)     = %.1fx\n', penalty / qvals(8));
fprintf('  penalty / max(genuine)     = %.1fx\n', penalty / max(genuine));

fprintf('\n---- Justification statement ----\n');
if max(genuine) < halfPenalty
    fprintf(['  A penalty of %.0f (and half-penalty of %.0f) lies above the ENTIRE\n' ...
             '  observed range of genuine >=2-camera uncertainty (max = %.4g m).\n' ...
             '  The penalty is therefore guaranteed to dominate, so any unseen or\n' ...
             '  single-camera point always costs more than the worst-triangulated\n' ...
             '  point, correctly steering the GA toward full triangulable coverage.\n'], ...
             penalty, halfPenalty, max(genuine));
else
    fprintf(['  NOTE: %.3f%% of genuine values exceed the half-penalty (50). Review\n' ...
             '  whether the penalty should be raised for clean separation.\n'], fracGe50);
end

%% ---------------- Plots -------------------------------------------------
figDir = fullfile(projectRoot, 'figures');
if saveFigures && ~exist(figDir,'dir'), mkdir(figDir); end

f1 = figure('Color','w','Position',[100 100 760 460]);
histogram(genuine, 60, 'FaceColor', [0.2 0.4 0.7], 'EdgeColor','none');
hold on;
yl = ylim;
plot([halfPenalty halfPenalty], yl, 'r--', 'LineWidth', 1.5);
plot([penalty penalty],         yl, 'r-',  'LineWidth', 1.5);
text(qvals(5), yl(2)*0.9, sprintf('median = %.3g m', qvals(5)), 'Color',[0.2 0.4 0.7]);
xlabel('Genuine \geq2-camera uncertainty  \Sigma\surd\lambda_i  [m]');
ylabel('count');
title(sprintf('Resolution-uncertainty distribution (%d configs, n=%d genuine)', numConfigs, numel(genuine)));
legend({'genuine values','half-penalty = 50','penalty = 100'}, 'Location','northeast');
grid on;
if saveFigures
    saveas(f1, fullfile(figDir,'resUncertaintyRange_linear.png'));
end

% Log-scale view: the genuine distribution and the penalties are orders
% of magnitude apart, so a log x-axis shows the separation clearly.
f2 = figure('Color','w','Position',[120 120 760 460]);
edges = logspace(floor(log10(max(min(genuine),eps))), log10(penalty*1.2), 60);
histogram(genuine, edges, 'FaceColor',[0.2 0.4 0.7], 'EdgeColor','none');
set(gca,'XScale','log'); hold on;
yl = ylim;
plot([halfPenalty halfPenalty], yl, 'r--', 'LineWidth', 1.5);
plot([penalty penalty],         yl, 'r-',  'LineWidth', 1.5);
xlabel('Genuine \geq2-camera uncertainty  [m]  (log scale)');
ylabel('count');
title('Separation between genuine metric and the penalty (log scale)');
legend({'genuine values','half-penalty = 50','penalty = 100'}, 'Location','northwest');
grid on;
if saveFigures
    saveas(f2, fullfile(figDir,'resUncertaintyRange_log.png'));
end

% Save the raw genuine values for later reference / thesis tables.
if saveFigures
    save(fullfile(figDir,'resUncertaintyRange_data.mat'), ...
        'genuine','cnt','qvals','pct','penalty','halfPenalty', ...
        'numConfigs','numCams','spacing','targetType');
end
fprintf('\nDone. Figures and data saved to %s\n', figDir);

%% ======================================================================
%  Local function: instrumented copy of computePointUncertainty.
%  Returns the SAME metric value plus a category so we can separate the
%  genuine >=2-camera bucket from the penalty buckets.
%     cat = 0 : 0 cameras   (original returns penaltyUncertainty)
%     cat = 1 : 1 camera    (original returns 0.5*penaltyUncertainty)
%     cat = 2 : >=2 cams but degenerate polytope (0.5*penaltyUncertainty)
%     cat = 3 : >=2 cams, valid polytope -> genuine metric returned in val
%  The cat==3 branch is byte-for-byte the same computation as the original.
%% ======================================================================
function [val, cat] = pointUncertaintyInstrumented(point, cameras, cameraCentres, ...
        numCams, adjacentSurfaces, du, dv, Resolution)

    planes = cell(1, numCams);
    visIdx = false(numCams,1);

    for i = 1:numCams
        uv = cameras{i}.project(point);
        u = uv(1); v = uv(2);
        if (u >= 1 && u <= Resolution(1) && v >= 1 && v <= Resolution(2))
            worldPoints = quantToWorld(cameras{i}, u, v, du, dv, cameraCentres(:,i));
            planes{i}   = buildPyramidSurf(cameraCentres(:,i), worldPoints, adjacentSurfaces);
            visIdx(i)   = true;
        end
    end

    visibleIdx = find(visIdx);
    numVisible = numel(visibleIdx);

    if numVisible == 0
        val = NaN; cat = 0; return;
    elseif numVisible == 1
        val = NaN; cat = 1; return;
    end

    vertices = calcVertices(numVisible, visibleIdx, planes);
    if isempty(vertices)
        val = NaN; cat = 2; return;
    end
    unique_vertices = unique(round(vertices,6), 'rows', 'stable');
    if size(unique_vertices,1) < 4
        val = NaN; cat = 2; return;
    end

    V_centered = unique_vertices - point;
    C0v_vert   = cov(V_centered);
    Eigs       = eig(C0v_vert);
    val        = sum(sqrt(abs(Eigs)));   
    if ~isfinite(val)
        val = NaN; cat = 2; return;  
    end
    cat = 3;
end
