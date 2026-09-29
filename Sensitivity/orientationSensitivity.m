function sweep = orientationSensitivity(varargin)
%ORIENTATIONSENSITIVITY  Roll-canonicalisation and angular-quantisation study.
%
%   sweep = orientationSensitivity() runs both halves of the orientation
%   sensitivity analysis on the Optimised GA Rig CF3 configuration, at the UAV target
%   space with 1 m grid spacing.
%
%   ONLY THE GA CONFIGURATION IS SNAPPED. The Manually Posed Rig is
%   already mounted in the room — its poses are measurements of hardware
%   that exists, not a specification anyone has to realise, so asking what
%   it would cost to build it on a 15-degree grid is a question about
%   nothing. The optimised layout is the one that has to be physically set
%   up with real brackets, and it is therefore the only one whose
%   mounting precision matters. The ad-hoc rig's CF3 is still evaluated
%   once, unsnapped, and carried through as a benchmark: it is the bar the
%   snapped optimum has to stay under to be worth building.
%
%   It answers two questions the raw GA output leaves open.
%
%   PART A — "how many cameras came out upside down, and what does fixing
%   it cost?"  The GA never constrains roll about the optical axis, so
%   optimal solutions routinely hang cameras inverted. uprightCameras
%   corrects this by adding 180 degrees to the roll gene (gamma), which
%   spins the image without moving the optical axis at all. A 180-degree
%   image rotation maps the sensor rectangle onto itself, so the field of
%   view is unchanged and the cost should not move.
%
%   It moves by a hair, and the reason is worth knowing. The visibility
%   test accepts u in [1, W] and v in [1, H], whose centre is
%   ((W+1)/2, (H+1)/2) = (640.5, 512.5), but the principal point is set to
%   (W/2, H/2) = (640, 512). The accepted band is therefore half a pixel
%   off-centre, so a 180-degree roll shifts it by one pixel and a target
%   sitting in the outermost pixel row or column can change state. In the
%   7-camera UAV case this is one (point, camera) pair out of 2835 and it
%   moves CF1 by about 5e-6. This part measures that delta rather than
%   assuming it away, and also prices the stricter alternative of levelling
%   every camera to the horizon, which IS a real change because the sensor
%   is 1280 x 1024 rather than square.
%
%   PART B — "the GA reports angles like 1.5708 and 0.7846 rad; what
%   happens if they are snapped to angles a person can actually set?"
%   Orientations are quantised to a grid of candidate step sizes (15
%   degrees among them) and the configuration is re-costed at each. Three
%   gene sets are priced separately, because they are not equivalent:
%
%     all       alpha, beta and gamma snapped
%     pointing  alpha and beta only — these steer the optical axis, so
%               this is the variant that actually moves coverage
%     roll      gamma only — a pure roll about the optical axis, which
%               should be very nearly free and acts as the control
%
%   Both parts run on the SAME target space and hardware specs as
%   runCameraOptimiser, so every number here is directly comparable to a
%   logged BestCost.
%
%   MONOTONICITY. Do not expect the penalty to grow smoothly with the step.
%   Snapping is not a small smooth perturbation, and CF2 counts points
%   crossing discrete triangulation-angle thresholds, so the cost surface is
%   genuinely rough at these scales. A coarser grid landing below a finer
%   one is that roughness, not an error.
%
%   READING THE PERCENTAGES. CF3 is utopia-shifted, so a well-optimised
%   configuration sits very close to zero and small absolute changes look
%   enormous in relative terms. Both the absolute and the relative change
%   are reported, and the ad-hoc rig's CF3 is carried through as a
%   real-world yardstick: a snapped GA layout that is still far below the
%   ad-hoc rig is still a good layout, whatever the percentage says.
%
%   USAGE
%       sweep = orientationSensitivity()                       % UAV
%       sweep = orientationSensitivity('TargetType', 2)        % UGV
%       sweep = orientationSensitivity('Steps', [5 10 15 30])
%       sweep = orientationSensitivity('SnapVariants', {'all'})
%       sweep = orientationSensitivity('AdHocReference', false)
%       sweep = orientationSensitivity('Plot', false)
%
%   Name-Value parameters
%     'NumCameras'    camera count. Default 7.
%     'TargetType'    1 = UAV (full volume), 2 = UGV (floor slab). Default 1.
%     'GridMode'      1 = uniform, 2 = centre-concentrated. Default 1.
%     'Spacing'       x-y evaluation grid spacing in m. Default 1.0.
%     'UGVmaxHeight'  UGV slab height in m. Default 0.5.
%     'UGVzSpacing'   UGV slab z step in m. Default 0.25.
%     'WeightUnc'     CF3 weight on resolution uncertainty. Default 0.5.
%     'WeightOcc'     CF3 weight on dynamic occlusion. Default 0.5.
%     'Steps'         orientation grids to test, in degrees.
%                     Default [1 2 5 10 15 22.5 30 45].
%     'SnapVariants'  any of {'all','pointing','roll'}. Default all three.
%     'PositionStep'  position grid in m applied alongside each orientation
%                     grid. Default [] — orientation only, which isolates
%                     the effect the analysis is about.
%     'TolerancePct'  CF3 deviation (%) above which a grid is called too
%                     coarse. Default 2.
%     'AdHocReference'  evaluate the Manually Posed Rig once (unsnapped)
%                     as a benchmark for the printed summary. It is never
%                     snapped and never plotted as a series — see above.
%                     Default true.
%     'Plot'          produce and save figures. Default true.
%
%   OUTPUT / SIDE EFFECTS
%   Returns the sweep struct, also written to
%     Results/Sensitivity/Orientation/orientation_sweep_<mode>_<ts>.mat
%   and, with 'Plot' true, four figures under
%     figures/Sensitivity/Orientation/
%
%   See also: uprightCameras, snapChromosome, cameraOrientationInfo,
%             spacingSensitivity_UAV, plotOrientationSensitivity.

    p = inputParser;
    addParameter(p, 'NumCameras',       7,     @isnumeric);
    addParameter(p, 'TargetType',       1,     @isnumeric);
    addParameter(p, 'GridMode',         1,     @isnumeric);
    addParameter(p, 'Spacing',          1.0,   @isnumeric);
    addParameter(p, 'UGVmaxHeight',     0.5,   @isnumeric);
    addParameter(p, 'UGVzSpacing',      0.25,  @isnumeric);
    addParameter(p, 'WeightUnc',        0.5,   @isnumeric);
    addParameter(p, 'WeightOcc',        0.5,   @isnumeric);
    addParameter(p, 'Steps',            [1 2 5 10 15 22.5 30 45], @isnumeric);
    addParameter(p, 'SnapVariants',     {'all','pointing','roll'}, @iscellstr);
    addParameter(p, 'PositionStep',     [],    @isnumeric);
    addParameter(p, 'TolerancePct',     2,     @isnumeric);
    addParameter(p, 'AdHocReference',   true,  @islogical);
    addParameter(p, 'Plot',             true,  @islogical);
    addParameter(p, 'Preset',           'optitrack_lab', @ischar);   % see runConfig
    parse(p, varargin{:});
    opts = p.Results;

    projectRoot = addProjectPaths();
    numCams = opts.NumCameras;

    if opts.TargetType == 2
        modeTag = 'UGV';
    else
        modeTag = 'UAV';
    end

    %% Snap variants -----------------------------------------------------
    allVariants = struct( ...
        'key',   {'all',              'pointing',                 'roll'}, ...
        'genes', {[1 2 3],            [1 2],                      3}, ...
        'label', {'all three angles', 'pointing only (alpha,beta)','roll only (gamma)'});
    keep = ismember({allVariants.key}, lower(opts.SnapVariants));
    assert(any(keep), 'orientationSensitivity:BadVariant', ...
        'SnapVariants must be a subset of {''all'',''pointing'',''roll''}.');
    variants = allVariants(keep);
    nV = numel(variants);

    %% Specs and search bounds — same builder as the GA runs, so snapping
    %% can never leave the feasible set
    run = runConfig(opts.Preset, 'UGV_MaxHeight', opts.UGVmaxHeight, ...
        'UGV_ZSpacing', opts.UGVzSpacing, 'Weights', [opts.WeightUnc, opts.WeightOcc]);
    [specs, problem] = buildRunSpecs(run, numCams, 3, opts.TargetType, opts.GridMode, opts.Spacing);
    volume = run.Volume;
    if opts.TargetType == 2
        volume(3, :) = [0, run.UGV_MaxHeight];
    end
    VarMin = problem.VarMin;
    VarMax = problem.VarMax;

    %% Configurations ---------------------------------------------------
    [gaChrom, bestRun] = loadBestCF3Config(numCams, opts.TargetType, opts.GridMode);
    configs(1).name       = sprintf('Optimised GA Rig CF3 (logged cost=%.5f)', bestRun.BestCost);
    configs(1).shortName  = 'Optimised GA Rig CF3';
    configs(1).chromosome = gaChrom;
    nC = numel(configs);

    % Benchmark only: the ad-hoc rig is already installed, so it is
    % evaluated once as it stands and never snapped.
    adhocCF3 = NaN;
    if opts.AdHocReference && numCams == 7
        adhocCost = evalChrom(buildOptiTrackChromosome(), specs, numCams);
        adhocCF3  = adhocCost.CF3;
    end

    steps = sort(opts.Steps(:).');
    nS    = numel(steps);

    %% Banner ------------------------------------------------------------
    fprintf('\n============================================================\n');
    fprintf('  Orientation sensitivity — %s, %d cams\n', modeTag, numCams);
    fprintf('============================================================\n');
    fprintf('Volume   : x=[%g %g], y=[%g %g], z=[%g %g] m\n', volume.');
    fprintf('Grid     : mode %d, spacing %g m -> %d target points\n', ...
        opts.GridMode, opts.Spacing, specs.NumPoints);
    fprintf('Weights  : w_unc=%.2f, w_occ=%.2f   (norm source: %s)\n', ...
        opts.WeightUnc, opts.WeightOcc, specs.PreComputed.normSource);
    fprintf('Steps    : '); fprintf('%g ', steps); fprintf('deg\n');
    fprintf('Variants : '); fprintf('%s ', variants.key); fprintf('\n');
    if ~isempty(opts.PositionStep)
        fprintf('Position : snapped to %g m alongside each orientation grid\n', opts.PositionStep);
    end
    fprintf('Config   : %s\n', configs(1).name);
    if isfinite(adhocCF3)
        fprintf('Benchmark: Manually Posed Rig, CF3 = %.5f (as installed, never snapped)\n', ...
            adhocCF3);
    end
    fprintf('------------------------------------------------------------\n');

    sweepStart = tic;

    blankUp = struct('config','', 'variant','', 'cost',[], 'report',[], ...
                     'info',[], 'chromosome',[]);
    upright = repmat(blankUp, nC*3, 1);
    uIdx = 0;

    blankSnap = struct('config','', 'variant','', 'stepDeg',NaN, 'cost',[], ...
                       'snapReport',[], 'chromosome',[]);
    snapRows = repmat(blankSnap, nC*nS*nV, 1);
    sIdx = 0;

    for c = 1:nC
        cname = configs(c).shortName;
        base  = configs(c).chromosome(:).';

        fprintf('\n=== Config %d/%d: %s ===\n', c, nC, configs(c).name);

        %% -------- PART A: roll canonicalisation ------------------------
        baseInfo = cameraOrientationInfo(base, numCams);
        nInv = nnz(baseInfo.Inverted & ~baseInfo.Degenerate);

        fprintf('\n-- Part A: roll about the optical axis --\n');
        fprintf('  %d of %d cameras are mounted upside down.\n', nInv, numCams);
        fprintf('  %-5s  %-12s  %-12s  %-10s\n', 'Cam', 'Tilt [deg]', 'Elev [deg]', 'Status');
        fprintf('  %s\n', repmat('-', 1, 48));
        for k = 1:numCams
            if baseInfo.Degenerate(k)
                status = 'vertical axis';
            elseif baseInfo.Inverted(k)
                status = 'INVERTED';
            else
                status = 'upright';
            end
            fprintf('  %-5d  %12.2f  %12.2f  %-10s\n', ...
                k, baseInfo.TiltDeg(k), baseInfo.ElevationDeg(k), status);
        end

        [flipChrom,  flipRep]  = uprightCameras(base, 'Mode', 'flip',  'NumCams', numCams);
        [levelChrom, levelRep] = uprightCameras(base, 'Mode', 'level', 'NumCams', numCams);

        baseCost  = evalChrom(base,       specs, numCams);
        flipCost  = evalChrom(flipChrom,  specs, numCams);
        levelCost = evalChrom(levelChrom, specs, numCams);

        vNames   = {'baseline', 'upright-flip', 'upright-level'};
        vCosts   = {baseCost, flipCost, levelCost};
        vReports = {[], flipRep, levelRep};
        vChroms  = {base, flipChrom, levelChrom};
        for v = 1:3
            uIdx = uIdx + 1;
            upright(uIdx).config     = cname;
            upright(uIdx).variant    = vNames{v};
            upright(uIdx).cost       = vCosts{v};
            upright(uIdx).report     = vReports{v};
            upright(uIdx).info       = cameraOrientationInfo(vChroms{v}, numCams);
            upright(uIdx).chromosome = vChroms{v};
        end

        fprintf('\n  %-16s %-11s %-11s %-11s %-12s %-9s\n', ...
            'Variant', 'CF1', 'CF2', 'CF3', 'dCF3', '2+ cams %');
        fprintf('  %s\n', repmat('-', 1, 74));
        for v = 1:3
            printCostRow(vNames{v}, vCosts{v}, baseCost);
        end
        fprintf(['\n  Flip changes gamma by 180 deg on %d camera(s) and nothing else, so\n' ...
                 '  the optical axis and the field of view are identical. The residual\n' ...
                 '  dCF3 above is the half-pixel asymmetry described in the header, not\n' ...
                 '  a change of coverage. Levelling rotates the non-square sensor and is\n' ...
                 '  the variant that can genuinely move the cost.\n'], flipRep.NumChanged);

        %% -------- PART B: angular quantisation -------------------------
        % Snap the UPRIGHT (flipped) configuration: that is the one a
        % technician would actually mount, so the quantisation penalty is
        % measured on top of the fix, not instead of it.
        fprintf('\n-- Part B: snapping orientations to a mounting grid --\n');
        fprintf('  Reference (unsnapped, upright): CF3 = %.5f\n', flipCost.CF3);

        for v = 1:nV
            fprintf('\n  [%s] %s\n', variants(v).key, variants(v).label);
            fprintf('  %-8s %-10s %-10s %-10s %-10s %-9s %-8s %-8s %-9s\n', ...
                'Step', 'CF1', 'CF2', 'CF3', 'dCF3', 'dCF3 %', 'maxAxis', 'maxRoll', '2+ cams %');
            fprintf('  %s\n', repmat('-', 1, 96));

            for s = 1:nS
                [snapped, snapRep] = snapChromosome(flipChrom, ...
                    'OrientationStepDeg', steps(s), ...
                    'Genes',              variants(v).genes, ...
                    'PositionStep',       opts.PositionStep, ...
                    'NumCams',            numCams, ...
                    'VarMin', VarMin, 'VarMax', VarMax);

                cost = evalChrom(snapped, specs, numCams);

                sIdx = sIdx + 1;
                snapRows(sIdx).config     = cname;
                snapRows(sIdx).variant    = variants(v).key;
                snapRows(sIdx).stepDeg    = steps(s);
                snapRows(sIdx).cost       = cost;
                snapRows(sIdx).snapReport = snapRep;
                snapRows(sIdx).chromosome = snapped;

                dAbs = cost.CF3 - flipCost.CF3;
                dPct = 100 * dAbs / abs(flipCost.CF3);
                fprintf('  %6.1f   %9.5f  %9.5f  %9.5f  %+9.5f  %+8.2f  %7.2f  %7.2f  %8.1f\n', ...
                    steps(s), cost.CF1, cost.CF2, cost.CF3, dAbs, dPct, ...
                    max(snapRep.AxisShiftDeg), ...
                    max(abs(snapRep.RollShiftDeg), [], 'omitnan'), cost.TwoPlusPct);
            end
        end

        % How close was the unconstrained optimum to the 15-degree grid
        % already? This is the "a lot of them were near pi, pi/2 and pi/4"
        % observation, quantified.
        printGridProximity(flipChrom, numCams, 15);
    end

    upright  = upright(1:uIdx);
    snapRows = snapRows(1:sIdx);
    sweepElapsed = toc(sweepStart);
    fprintf('\nOrientation sweep complete in %.1f s.\n', sweepElapsed);

    %% Pack --------------------------------------------------------------
    sweep = struct();
    sweep.modeTag       = modeTag;
    sweep.timestamp     = datetime('now');
    sweep.numCams       = numCams;
    sweep.targetType    = opts.TargetType;
    sweep.gridMode      = opts.GridMode;
    sweep.spacing       = opts.Spacing;
    sweep.volume        = volume;
    sweep.numPoints     = specs.NumPoints;
    sweep.weightUnc     = opts.WeightUnc;
    sweep.weightOcc     = opts.WeightOcc;
    sweep.steps         = steps;
    sweep.variants      = variants;
    sweep.positionStep  = opts.PositionStep;
    sweep.tolerancePct  = opts.TolerancePct;
    sweep.configs       = configs;
    sweep.adhocCF3      = adhocCF3;
    sweep.upright       = upright;
    sweep.snap          = snapRows;
    sweep.uncertNorm    = specs.PreComputed.uncertNorm;
    sweep.occlNorm      = specs.PreComputed.occlNorm;
    sweep.normSource    = specs.PreComputed.normSource;
    sweep.elapsed       = sweepElapsed;

    %% Save --------------------------------------------------------------
    outDir = fullfile(projectRoot, 'Results', 'Sensitivity', 'Orientation');
    if ~isfolder(outDir), mkdir(outDir); end
    ts = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    matFile = fullfile(outDir, sprintf('orientation_sweep_%s_%s.mat', modeTag, ts));
    save(matFile, 'sweep');
    fprintf('Sweep saved to: %s\n', matFile);

    %% Recommendation ----------------------------------------------------
    printRecommendation(sweep);

    %% Plot --------------------------------------------------------------
    if opts.Plot
        plotOrientationSensitivity(sweep, projectRoot, ts);
    end
end


%% ======================================================================
%  LOCAL HELPERS
%  ======================================================================

function cost = evalChrom(chrom, specs, numCams)
% Evaluate one chromosome exactly the way the GA's cost function does,
% plus a coverage summary so the "roll does not change coverage" claim is
% checked rather than asserted.
    [cameras, camCentres] = setupCameras(chrom, numCams, specs.Resolution, ...
        specs.Focal, specs.FocalWide, specs.PrincipalPoint, specs.PixelSize);

    cost.CF1 = resUncertainty(specs, cameras, camCentres);
    cost.CF2 = dynamicOcclusion(specs, cameras, camCentres);
    [cost.CF3, cost.Junc, cost.Jocc] = cf3Terms(cost.CF1, cost.CF2, specs);

    [~, covStats] = perTargetCoverage(chrom, specs);
    cost.MeanCams   = covStats.avg;
    cost.ZeroPct    = covStats.zeroPct;
    cost.TwoPlusPct = covStats.twoPlusPct;
end


function printCostRow(label, cost, refCost)
    fprintf('  %-16s %10.5f  %10.5f  %10.5f  %+11.3e  %8.1f\n', ...
        label, cost.CF1, cost.CF2, cost.CF3, cost.CF3 - refCost.CF3, cost.TwoPlusPct);
end


function printGridProximity(chrom, numCams, stepDeg)
% Distance from every orientation gene to the nearest multiple of stepDeg,
% before any snapping. Small residuals mean the GA had already converged
% onto near-mountable angles on its own.
    res   = zeros(numCams, 3);
    names = {'alpha','beta','gamma'};
    for c = 1:numCams
        idx = (c-1)*6 + 4;
        ang = rad2deg(chrom(idx:idx+2));
        res(c,:) = abs(ang - round(ang / stepDeg) * stepDeg);
    end
    allRes = res(:);
    fprintf('\n  Proximity of the unsnapped optimum to the %g deg grid:\n', stepDeg);
    fprintf('    mean %.2f deg, median %.2f deg, max %.2f deg  (worst case is %.1f deg)\n', ...
        mean(allRes), median(allRes), max(allRes), stepDeg/2);
    fprintf('    %d of %d genes (%.0f%%) already lie within 2 deg of a grid angle.\n', ...
        nnz(allRes <= 2), numel(allRes), 100*nnz(allRes <= 2)/numel(allRes));
    for g = 1:3
        fprintf('    %-6s mean %.2f deg\n', names{g}, mean(res(:,g)));
    end
    fprintf(['    A uniformly distributed angle would average %.2f deg from the grid,\n' ...
             '    so compare against that before reading anything into these.\n'], stepDeg/4);
end


function printRecommendation(sweep)
% Coarsest grid whose CF3 penalty stays inside the tolerance.
    if isempty(sweep.tolerancePct), return; end

    fprintf('\n------------------------------------------------------------\n');
    fprintf('  Coarsest mountable grid within %.1f%% of the unsnapped CF3\n', sweep.tolerancePct);
    fprintf('------------------------------------------------------------\n');

    % Ad-hoc rig CF3 as a real-world yardstick for the absolute numbers.
    adhoc = NaN;
    if isfield(sweep, 'adhocCF3'), adhoc = sweep.adhocCF3; end

    cfgNames = unique({sweep.snap.config}, 'stable');
    for c = 1:numel(cfgNames)
        ref = refCF3(sweep, cfgNames{c});
        fprintf('\n  %s   (unsnapped CF3 = %.5f)\n', cfgNames{c}, ref);

        for v = 1:numel(sweep.variants)
            key  = sweep.variants(v).key;
            rows = sweep.snap(strcmp({sweep.snap.config}, cfgNames{c}) & ...
                              strcmp({sweep.snap.variant}, key));
            if isempty(rows), continue; end
            [~, order] = sort([rows.stepDeg]);
            rows = rows(order);

            cf3 = arrayfun(@(r) r.cost.CF3, rows);
            dev = 100 * (cf3 - ref) / abs(ref);
            ok  = abs(dev) <= sweep.tolerancePct;

            if any(ok)
                best = find(ok, 1, 'last');
                fprintf('    %-10s coarsest within tolerance: %g deg (CF3 %+.2f%%)\n', ...
                    key, rows(best).stepDeg, dev(best));
            else
                fprintf('    %-10s no tested grid stays within tolerance\n', key);
            end

            % The finest grid tested is the practical resolution of this
            % comparison. CF2 counts points crossing discrete
            % triangulation-angle thresholds, so the cost surface is
            % genuinely rough at small angular scales — a deviation of the
            % same size as this floor says nothing about the grid.
            fprintf('    %-10s floor: the finest grid tested (%g deg) already moves CF3 %+.2f%%\n', ...
                '', rows(1).stepDeg, dev(1));

            at15 = find([rows.stepDeg] == 15, 1);
            if ~isempty(at15)
                msg = sprintf(['    %-10s at 15 deg: CF3 %.5f (%+.2f%%), ' ...
                               'max axis shift %.1f deg, max roll shift %.1f deg'], ...
                    '', cf3(at15), dev(at15), ...
                    max(rows(at15).snapReport.AxisShiftDeg), ...
                    max(abs(rows(at15).snapReport.RollShiftDeg), [], 'omitnan'));
                fprintf('%s\n', msg);
                if isfinite(adhoc)
                    fprintf(['    %-10s              = %.2f%% of the installed ad-hoc ' ...
                             'rig''s CF3 (%.5f)\n'], '', 100*cf3(at15)/adhoc, adhoc);
                end
            end
        end
    end
    fprintf('\n');
end


function ref = refCF3(sweep, configName)
% The unsnapped, upright configuration is the reference for Part B.
    rows = sweep.upright(strcmp({sweep.upright.config}, configName) & ...
                         strcmp({sweep.upright.variant}, 'upright-flip'));
    ref = rows(1).cost.CF3;
end
