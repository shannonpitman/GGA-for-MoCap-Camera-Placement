function S = mountingSensitivity(varargin)
%MOUNTINGSENSITIVITY  How much does real mounting error cost an optimised rig?
%
%   S = mountingSensitivity()                     % best CF3 rig, UAV uniform
%   S = mountingSensitivity('TargetType', 2)      % UGV slab
%   S = mountingSensitivity('Quick', true)        % smoke test
%
%   The optimised rig is uprighted, then perturbed the way installation
%   actually perturbs it, and re-costed with the GA's own cost functions.
%
%   PART A - one gene at a time. Each camera's genes are perturbed alone by
%   +/- each level, and the mean |change in cost| over cameras and signs is
%   recorded. Genes:
%     x, y, z         world position [cm]; the camera is then projected back
%                     onto its mount region (a wall camera cannot move into
%                     the wall), as in the GA
%     pan, tilt, roll physical pan-tilt-head angles [deg]: pan about world
%                     vertical, tilt about the camera's image-right axis,
%                     roll about its optical axis
%     rx, ry, rz      raw rotation-vector genes [deg], for reference
%   Genes are ranked by their effect at the achievable tolerance.
%
%   PART B - Monte Carlo. Every camera is perturbed at once with Gaussian
%   error, N draws per level, for three groups: position only (x,y,z),
%   orientation only (pan,tilt,roll) and both. Reports median and 90th
%   percentile of J, J_res, J_occ and the % of points seen by >= 2 cameras.
%
%   PART C - hypothesis test at the achievable tolerance
%   (AchievablePos, AchievableAng; default 1 cm and 1 deg = tape measure +
%   pan-tilt head, about +/-2 cm / +/-2 deg). Draws are paired (same seed per
%   draw), and the test statistic is d = dJ(orientation) - dJ(position).
%   Reported: median d with a bootstrap 95% CI, share of draws with d > 0,
%   and a two-sided sign-test p-value. "Orientation matters more" is
%   supported when the CI lies above zero.
%
%   Results go to Results/Sensitivity/Mounting/ (.mat + .txt) and figures
%   to figures/Sensitivity/Mounting/.
%
%   Name-value parameters
%     'Chromosome'     rig to test (default: best logged CF3 rig)
%     'Preset'         runConfig preset                 'optitrack_lab'
%     'NumCameras'     7   'TargetType' 1   'GridMode' 1   'Spacing' 1.0
%     'PosLevelsCm'    [0.5 1 2 5 10]
%     'AngLevelsDeg'   [0.5 1 2 5 10]
%     'AchievablePos'  1   [cm, 1 s.d.]
%     'AchievableAng'  1   [deg, 1 s.d.]
%     'NumDraws'       200
%     'Seed'           1
%     'Plot'           true
%     'Quick'          false

    projectRoot = addProjectPaths();

    p = inputParser;
    addParameter(p, 'Chromosome',    [],   @isnumeric);
    addParameter(p, 'Preset',        'optitrack_lab', @ischar);
    addParameter(p, 'NumCameras',    7,    @isnumeric);
    addParameter(p, 'TargetType',    1,    @isnumeric);
    addParameter(p, 'GridMode',      1,    @isnumeric);
    addParameter(p, 'Spacing',       1.0,  @isnumeric);
    addParameter(p, 'PosLevelsCm',   [0.5 1 2 5 10], @isnumeric);
    addParameter(p, 'AngLevelsDeg',  [0.5 1 2 5 10], @isnumeric);
    addParameter(p, 'AchievablePos', 1,    @isnumeric);
    addParameter(p, 'AchievableAng', 1,    @isnumeric);
    addParameter(p, 'NumDraws',      200,  @isnumeric);
    addParameter(p, 'Seed',          1,    @isnumeric);
    addParameter(p, 'Plot',          true, @islogical);
    addParameter(p, 'Quick',         false, @islogical);
    addParameter(p, 'OutputDir',     '',   @ischar);
    parse(p, varargin{:});
    o = p.Results;

    if o.Quick
        o.Spacing = 2.0;  o.NumDraws = 6;  o.PosLevelsCm = [1 5];  o.AngLevelsDeg = [1 5];
    end
    if isempty(o.OutputDir)
        o.OutputDir = fullfile(projectRoot, 'Results', 'Sensitivity', 'Mounting');
    end
    if ~isfolder(o.OutputDir), mkdir(o.OutputDir); end

    %% Setup: same specs as the GA runs
    cfg = runConfig(o.Preset);
    specs = buildRunSpecs(cfg, o.NumCameras, 3, o.TargetType, o.GridMode, o.Spacing);
    if isempty(o.Chromosome)
        o.Chromosome = loadBestCF3Config(o.NumCameras, o.TargetType, o.GridMode);
    end
    base = uprightCameras(enforceMount(o.Chromosome, specs.MountRegions));
    nCam = o.NumCameras;

    ref = evaluateRig(base, specs);
    fprintf('\n  mountingSensitivity — %s, TT%d GM%d, %d points\n', o.Preset, ...
        o.TargetType, o.GridMode, specs.NumPoints);
    fprintf('  Unperturbed: J = %.5f | J_res = %.5f m | J_occ = %.2f deg | >=2 cams %.1f%%\n\n', ...
        ref.J, ref.Jres, ref.Jocc, ref.TwoPlus);

    %% PART A: one gene at a time
    genes = {'x','y','z','pan','tilt','roll','rx','ry','rz'};
    isPos = [true true true false false false false false false];
    nG = numel(genes);
    oat = struct('Gene', genes, 'Levels', [], 'MeanAbsDJ', [], 'MeanAbsDJres', [], 'MeanAbsDJocc', []);

    for g = 1:nG
        if isPos(g), levels = o.PosLevelsCm; else, levels = o.AngLevelsDeg; end
        dJ = zeros(numel(levels), 1);  dR = dJ;  dO = dJ;
        for L = 1:numel(levels)
            jobs = zeros(0, 2);                          % [camera, sign]
            for c = 1:nCam, jobs = [jobs; c 1; c -1]; end %#ok<AGROW>
            vals = zeros(size(jobs, 1), 3);
            parfor k = 1:size(jobs, 1)
                delta = zeros(nCam, 9);
                delta(jobs(k,1), g) = jobs(k,2) * levels(L);
                e = evaluateRig(perturbRig(base, delta, specs.MountRegions), specs);
                vals(k, :) = [abs(e.J - ref.J), abs(e.Jres - ref.Jres), abs(e.Jocc - ref.Jocc)];
            end
            dJ(L) = mean(vals(:,1));  dR(L) = mean(vals(:,2));  dO(L) = mean(vals(:,3));
        end
        oat(g).Levels = levels;
        oat(g).MeanAbsDJ = dJ;
        oat(g).MeanAbsDJres = dR;
        oat(g).MeanAbsDJocc = dO;
        fprintf('  OAT %-5s ', genes{g});
        fprintf('%8.5f', dJ); fprintf('   (mean |dJ| per level)\n');
    end

    % Rank at the achievable tolerance
    achJ = zeros(nG, 1);
    for g = 1:nG
        if isPos(g), a = o.AchievablePos; else, a = o.AchievableAng; end
        achJ(g) = interp1(oat(g).Levels, oat(g).MeanAbsDJ, a, 'linear', 'extrap');
    end
    [~, rankIdx] = sort(achJ, 'descend');
    fprintf('\n  Ranking at achievable tolerance (%.3g cm, %.3g deg):\n', o.AchievablePos, o.AchievableAng);
    for k = rankIdx.'
        fprintf('    %-5s  mean |dJ| = %.5f\n', genes{k}, achJ(k));
    end

    %% PART B: Monte Carlo, all cameras at once
    groups = {'position', 'orientation', 'both'};
    nL = max(numel(o.PosLevelsCm), numel(o.AngLevelsDeg));
    posL = padLevels(o.PosLevelsCm, nL);
    angL = padLevels(o.AngLevelsDeg, nL);
    mc = struct([]);
    for gi = 1:numel(groups)
        for L = 1:nL
            [Js, Rs, Os, Ts] = monteCarlo(base, specs, groups{gi}, posL(L), angL(L), o.NumDraws, o.Seed);
            mc(end+1).Group = groups{gi}; %#ok<AGROW>
            mc(end).PosSigmaCm = posL(L) * ~strcmp(groups{gi}, 'orientation');
            mc(end).AngSigmaDeg = angL(L) * ~strcmp(groups{gi}, 'position');
            mc(end).J = Js;  mc(end).Jres = Rs;  mc(end).Jocc = Os;  mc(end).TwoPlus = Ts;
            fprintf('  MC %-11s pos %5.2f cm  ang %5.2f deg | J median %.5f  P90 %.5f | >=2 cams median %.1f%%\n', ...
                groups{gi}, mc(end).PosSigmaCm, mc(end).AngSigmaDeg, ...
                median(Js), prctile(Js, 90), median(Ts));
        end
    end

    %% PART C: paired hypothesis test at the achievable tolerance
    [JposA, RposA, OposA] = monteCarlo(base, specs, 'position',    o.AchievablePos, o.AchievableAng, o.NumDraws, o.Seed);
    [JangA, RangA, OangA] = monteCarlo(base, specs, 'orientation', o.AchievablePos, o.AchievableAng, o.NumDraws, o.Seed);
    H.J    = pairedTest(JangA - ref.J,    JposA - ref.J,    o.Seed);
    H.Jres = pairedTest(RangA - ref.Jres, RposA - ref.Jres, o.Seed);
    H.Jocc = pairedTest(OangA - ref.Jocc, OposA - ref.Jocc, o.Seed);
    fprintf('\n  Hypothesis: orientation error costs more than position error\n');
    fprintf('  at %.3g cm / %.3g deg (1 s.d.), %d paired draws\n', o.AchievablePos, o.AchievableAng, o.NumDraws);
    for f = {'J', 'Jres', 'Jocc'}
        h = H.(f{1});
        fprintf('    %-4s median d = %+.5f  95%% CI [%+.5f, %+.5f]  P(d>0) = %.2f  sign-test p = %.3g  -> %s\n', ...
            f{1}, h.MedianD, h.CI(1), h.CI(2), h.ShareAbove, h.SignP, h.Verdict);
    end

    %% Save
    S.Options = o;
    S.Reference = ref;
    S.Manual = [];
    if o.NumCameras == 7
        S.Manual = evaluateRig(buildOptiTrackChromosome(), specs);
    end
    S.OAT = oat;
    S.RankAchievable = table(genes(rankIdx).', achJ(rankIdx), 'VariableNames', {'Gene', 'MeanAbsDJ'});
    S.MonteCarlo = mc;
    S.Hypothesis = H;

    stamp = string(datetime('now'), 'yyyyMMdd_HHmmss');
    tag = sprintf('TT%d_GM%d_%dC', o.TargetType, o.GridMode, o.NumCameras);
    outFile = fullfile(o.OutputDir, sprintf('mountingSensitivity_%s_%s.mat', tag, stamp));
    save(outFile, 'S');
    writeSummary(strrep(outFile, '.mat', '.txt'), S, genes, rankIdx, achJ);
    fprintf('\n  Saved %s\n', outFile);

    if o.Plot
        plotMounting(S, genes, isPos, achJ, rankIdx, fullfile(projectRoot, 'figures', 'Sensitivity', 'Mounting'), tag);
    end
end

%% ========================================================================
function e = evaluateRig(chrom, specs)
    [cams, cc] = setupCameras(chrom, specs.Cams, specs.Resolution, specs.Focal, ...
        specs.FocalWide, specs.PrincipalPoint, specs.PixelSize);
    e.Jres = resUncertainty(specs, cams, cc);
    e.Jocc = dynamicOcclusion(specs, cams, cc);
    e.J    = cf3Terms(e.Jres, e.Jocc, specs);
    vis = projectVisibilityOcclusion(cams, specs.Target, cc, specs.Resolution, ...
        specs.PreComputed.maxCameraRange, specs.PreComputed.maxCameraRangeWide, specs.FocalWide);
    e.TwoPlus = 100 * mean(sum(vis, 2) >= 2);
end

function chrom = perturbRig(chrom, delta, regions)
% delta: numCams x 9 = [dx dy dz (cm), pan tilt roll (deg), drx dry drz (deg)]
    for c = 1:size(delta, 1)
        s = (c-1)*6;
        d = delta(c, :);
        if any(d(1:3))
            chrom(s+1:s+3) = projectToMount(chrom(s+1:s+3) + d(1:3)/100, regions);
        end
        if any(d(4:6))
            R = genesToRotm(chrom(s+4:s+6));
            R = rotZd(d(4)) * R * rotXd(d(5)) * rotZd(d(6));   % pan (world z), tilt (cam x), roll (cam z)
            chrom(s+4:s+6) = rotmToGenes(R);
        end
        if any(d(7:9))
            chrom(s+4:s+6) = rotvecWrap(chrom(s+4:s+6) + deg2rad(d(7:9)));
        end
    end
end

function R = rotXd(a)
    R = [1 0 0; 0 cosd(a) -sind(a); 0 sind(a) cosd(a)];
end

function R = rotZd(a)
    R = [cosd(a) -sind(a) 0; sind(a) cosd(a) 0; 0 0 1];
end

function [Js, Rs, Os, Ts] = monteCarlo(base, specs, group, posCm, angDeg, N, seed)
    nCam = specs.Cams;
    Js = zeros(N, 1);  Rs = Js;  Os = Js;  Ts = Js;
    regions = specs.MountRegions;
    parfor k = 1:N
        stream = RandStream('twister', 'Seed', seed*100000 + k);   % paired across groups
        zPos = randn(stream, nCam, 3);
        zAng = randn(stream, nCam, 3);
        delta = zeros(nCam, 9);
        if ~strcmp(group, 'orientation'), delta(:, 1:3) = posCm * zPos; end
        if ~strcmp(group, 'position'),    delta(:, 4:6) = angDeg * zAng; end
        e = evaluateRig(perturbRig(base, delta, regions), specs);
        Js(k) = e.J;  Rs(k) = e.Jres;  Os(k) = e.Jocc;  Ts(k) = e.TwoPlus;
    end
end

function h = pairedTest(dAng, dPos, seed)
    d = dAng(:) - dPos(:);
    h.MedianD = median(d);
    rs = RandStream('twister', 'Seed', seed);
    B = 2000;
    boot = zeros(B, 1);
    for b = 1:B
        boot(b) = median(d(randi(rs, numel(d), numel(d), 1)));
    end
    h.CI = prctile(boot, [2.5 97.5]);
    nz = d(d ~= 0);
    h.ShareAbove = mean(d > 0);
    h.SignP = min(1, 2 * binocdf(min(sum(nz > 0), sum(nz < 0)), numel(nz), 0.5));
    if h.CI(1) > 0
        h.Verdict = 'orientation worse';
    elseif h.CI(2) < 0
        h.Verdict = 'position worse';
    else
        h.Verdict = 'no clear difference';
    end
end

function L = padLevels(v, n)
    L = v(min(1:n, numel(v)));
end

function writeSummary(file, S, genes, rankIdx, achJ)
    fid = fopen(file, 'w');
    c = onCleanup(@() fclose(fid));
    o = S.Options;
    fprintf(fid, 'mountingSensitivity summary\n===========================\n');
    fprintf(fid, 'Preset %s, TT%d GM%d, spacing %.2f m, %d draws\n', o.Preset, ...
        o.TargetType, o.GridMode, o.Spacing, o.NumDraws);
    r = S.Reference;
    fprintf(fid, 'Unperturbed: J %.5f | J_res %.5f m | J_occ %.2f deg | >=2 cams %.1f%%\n', ...
        r.J, r.Jres, r.Jocc, r.TwoPlus);
    if ~isempty(S.Manual)
        m = S.Manual;
        fprintf(fid, 'Manual rig:  J %.5f | J_res %.5f m | J_occ %.2f deg | >=2 cams %.1f%%\n', ...
            m.J, m.Jres, m.Jocc, m.TwoPlus);
    end
    fprintf(fid, '\nGene ranking at achievable tolerance (%.3g cm / %.3g deg):\n', ...
        o.AchievablePos, o.AchievableAng);
    for k = rankIdx.'
        fprintf(fid, '  %-5s mean |dJ| %.5f\n', genes{k}, achJ(k));
    end
    fprintf(fid, '\nHypothesis (d = dJ_orientation - dJ_position, paired draws):\n');
    for f = {'J', 'Jres', 'Jocc'}
        h = S.Hypothesis.(f{1});
        fprintf(fid, '  %-4s median d %+.5f  CI [%+.5f, %+.5f]  P(d>0) %.2f  sign p %.3g  -> %s\n', ...
            f{1}, h.MedianD, h.CI(1), h.CI(2), h.ShareAbove, h.SignP, h.Verdict);
    end
    fprintf(fid, '\nMonte Carlo:\n');
    for k = 1:numel(S.MonteCarlo)
        m = S.MonteCarlo(k);
        fprintf(fid, '  %-11s pos %5.2f cm ang %5.2f deg | J median %.5f P90 %.5f | J_res median %.5f | J_occ median %.2f | >=2 cams %.1f%%\n', ...
            m.Group, m.PosSigmaCm, m.AngSigmaDeg, median(m.J), prctile(m.J, 90), ...
            median(m.Jres), median(m.Jocc), median(m.TwoPlus));
    end
end

function plotMounting(S, genes, isPos, achJ, rankIdx, figDir, tag)
    if ~isfolder(figDir), mkdir(figDir); end
    o = S.Options;
    fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 12 4.2]);
    tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact');

    ax = nexttile(tl); hold(ax, 'on');
    for g = find(isPos)
        plot(ax, S.OAT(g).Levels, S.OAT(g).MeanAbsDJ, '-o', 'LineWidth', 1.4, 'DisplayName', genes{g});
    end
    xline(ax, o.AchievablePos, '--k', 'achievable', 'HandleVisibility', 'off');
    set(ax, 'XScale', 'log'); xlabel(ax, 'Position error (cm)'); ylabel(ax, 'Mean |\DeltaJ|');
    title(ax, 'Position genes'); legend(ax, 'Location', 'southeast'); grid(ax, 'on');

    ax = nexttile(tl); hold(ax, 'on');
    for g = find(~isPos)
        plot(ax, S.OAT(g).Levels, S.OAT(g).MeanAbsDJ, '-o', 'LineWidth', 1.4, 'DisplayName', genes{g});
    end
    xline(ax, o.AchievableAng, '--k', 'achievable', 'HandleVisibility', 'off');
    set(ax, 'XScale', 'log'); xlabel(ax, 'Orientation error (deg)');
    title(ax, 'Orientation genes'); legend(ax, 'Location', 'southeast'); grid(ax, 'on');

    ax = nexttile(tl);
    bar(ax, categorical(genes(rankIdx), genes(rankIdx)), achJ(rankIdx));
    ylabel(ax, 'Mean |\DeltaJ| at achievable tolerance');
    title(ax, sprintf('Ranking (%.3g cm / %.3g deg)', o.AchievablePos, o.AchievableAng)); grid(ax, 'on');

    sty = gaPlotStyle();
    applyThesisStyle(fig);
    exportThesisFigure(fig, fullfile(figDir, sprintf('MountingSensitivity_%s', tag)), ...
        'Background', sty.ExportBgColor, 'Quiet', true);
end
