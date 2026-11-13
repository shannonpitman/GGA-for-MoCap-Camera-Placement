function T = reportABWarmStart(resultFile)
%REPORTABWARMSTART  Statistics and figures for an abWarmStart study.
%
%   T = reportABWarmStart('Results/Tuning/ABWarmStart/<stamp>/abWarmStart.mat')
%
%   Per camera count, for each arm (method x seed count), using completed
%   repeats only:
%     - best cost: median, IQR, mean, s.d., min
%     - share of runs that beat the best seed in the frozen archive
%     - layout spread: mean layoutDistance of each run's best rig to the
%       arm's best rig [m] (how different the warm solutions are)
%     - generations to come within 1% of the run's final cost (median)
%   Tests:
%     - Friedman test across all arms (repeats as blocks, complete repeats)
%     - paired Wilcoxon signed-rank: old vs new at each seed count, and
%       seed counts within each method
%     - Brown-Forsythe test of equal spread of cost across arms (a Levene
%       test on deviations from each arm's median; robust to non-normal
%       costs). Small p = the arms differ in how much the runs vary.
%   Writes report.txt, costs.png/.fig and convergence.png/.fig next to the
%   result file and returns the summary table.

    S = load(resultFile, 'R');
    R = S.R;
    outDir = fileparts(resultFile);
    runs = R.runs(~isnan([R.runs.BestCost]));
    if isempty(runs)
        error('reportABWarmStart:NoRuns', 'No completed runs in %s.', resultFile);
    end
    armName = @(r) sprintf('%s_%d', r.Method, r.NumSeeds);
    arms = unique(arrayfun(armName, runs, 'UniformOutput', false), 'stable');

    fid = fopen(fullfile(outDir, 'report.txt'), 'w');
    cleanup = onCleanup(@() fclose(fid));
    out = @(varargin) fprintfBoth(fid, varargin{:});
    out('Warm-start A/B study  %s\n', resultFile);
    out('Methods: %s | seed counts: %s | min seed distance %.2f m\n\n', ...
        strjoin(R.opts.Methods, ', '), mat2str(R.opts.SeedCounts), R.opts.MinDist);

    T = table();
    figC = figure('Name', 'Warm-start A/B: costs', 'Position', [100 100 560*numel(R.opts.Cameras) 420]);
    figV = figure('Name', 'Warm-start A/B: convergence', 'Position', [100 600 560*numel(R.opts.Cameras) 420]);
    for ci = 1:numel(R.opts.Cameras)
        nc = R.opts.Cameras(ci);
        rc = runs([runs.NumCams] == nc);
        if isempty(rc), continue; end
        snap = R.snapshot([R.snapshot.NumCams] == nc);
        bestSeed = snap.Pool(1).Cost;
        out('=== %d cameras  (%s)  archive: %d solutions, best %.5f ===\n', ...
            nc, snap.Key, numel(snap.Pool), bestSeed);

        % repeats x arms matrix (NaN where missing)
        nRep = max([rc.Repeat]);
        M = nan(nRep, numel(arms));
        for a = 1:numel(arms)
            ra = rc(strcmp(arrayfun(armName, rc, 'UniformOutput', false), arms{a}));
            M([ra.Repeat], a) = [ra.BestCost];

            costs = [ra.BestCost];
            C = vertcat(ra.BestChromosome);
            [~, ib] = min(costs);
            spread = arrayfun(@(k) layoutDistance(C(k, :), C(ib, :)), 1:numel(costs));
            spread(ib) = [];
            g1 = arrayfun(@(r) find(r.Convergence <= r.Convergence(end) * 1.01, 1), ra);

            row = table(nc, string(arms{a}), numel(costs), median(costs), iqr(costs), ...
                mean(costs), std(costs), min(costs), mean(costs < bestSeed - 1e-9), ...
                mean(spread), median(g1), ...
                'VariableNames', {'Cams', 'Arm', 'N', 'Median', 'IQR', 'Mean', 'SD', ...
                'Min', 'BeatBestSeed', 'LayoutSpread_m', 'GenTo1pct'});
            T = [T; row]; %#ok<AGROW>
        end
        out('%s\n', evalc('disp(T(T.Cams == nc, :))'));

        % Friedman over complete repeats
        full = all(~isnan(M), 2);
        if sum(full) >= 2 && numel(arms) > 1
            pF = friedman(M(full, :), 1, 'off');
            out('Friedman (%d complete repeats, %d arms): p = %.4g\n', sum(full), numel(arms), pF);
        end

        % Paired Wilcoxon: old vs new per seed count, seed counts per method
        for s = R.opts.SeedCounts
            pairTest(out, M, arms, sprintf('best_mutant_%d', s), sprintf('top_n_%d', s));
        end
        for m = R.opts.Methods
            sc = R.opts.SeedCounts;
            for k = 2:numel(sc)
                pairTest(out, M, arms, sprintf('%s_%d', m{1}, sc(1)), sprintf('%s_%d', m{1}, sc(k)));
            end
        end

        % Brown-Forsythe: equal spread of cost across arms
        y = M(:);  g = repelem(arms(:), nRep, 1);
        ok = ~isnan(y);
        if numel(unique(g(ok))) > 1
            pBF = vartestn(y(ok), g(ok), 'TestType', 'BrownForsythe', 'Display', 'off');
            out('Brown-Forsythe (equal spread of cost across arms): p = %.4g\n\n', pBF);
        end

        % Figures
        figure(figC); subplot(1, numel(R.opts.Cameras), ci);
        boxchart(categorical(g(ok), arms, 'Ordinal', true), y(ok)); hold on;
        yline(bestSeed, 'k--', 'best seed');
        title(sprintf('%d cameras', nc)); ylabel('final best cost (CF3)'); grid on;

        figure(figV); subplot(1, numel(R.opts.Cameras), ci); hold on;
        for a = 1:numel(arms)
            ra = rc(strcmp(arrayfun(armName, rc, 'UniformOutput', false), arms{a}));
            Cv = vertcat(ra.Convergence);
            plot(median(Cv, 1), 'LineWidth', 1.5, 'DisplayName', strrep(arms{a}, '_', ' '));
        end
        yline(bestSeed, 'k--', 'HandleVisibility', 'off');
        xlabel('generation'); ylabel('median best cost'); title(sprintf('%d cameras', nc));
        legend('Location', 'northeast'); grid on;
    end
    saveas(figC, fullfile(outDir, 'costs.png'));        savefig(figC, fullfile(outDir, 'costs.fig'));
    saveas(figV, fullfile(outDir, 'convergence.png'));  savefig(figV, fullfile(outDir, 'convergence.fig'));
    save(fullfile(outDir, 'summary.mat'), 'T');
    fprintf('Report written to %s\n', fullfile(outDir, 'report.txt'));
end

function pairTest(out, M, arms, a, b)
    ia = find(strcmp(arms, a));  ib = find(strcmp(arms, b));
    if isempty(ia) || isempty(ib), return; end
    ok = ~isnan(M(:, ia)) & ~isnan(M(:, ib));
    if sum(ok) < 2, return; end
    d = M(ok, ib) - M(ok, ia);
    pW = signrank(M(ok, ia), M(ok, ib));
    out('Wilcoxon %s vs %s (n=%d): median diff %+.5f, %s better in %d/%d, p = %.4g\n', ...
        b, a, sum(ok), median(d), b, sum(d < 0), sum(ok), pW);
end

function fprintfBoth(fid, varargin)
    fprintf(varargin{:});
    fprintf(fid, varargin{:});
end
