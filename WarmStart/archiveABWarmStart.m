function n = archiveABWarmStart(resultFile)
%ARCHIVEABWARMSTART  Add an abWarmStart study's final rigs to the seed archive.
%
%   n = archiveABWarmStart('Results/Tuning/ABWarmStart/<stamp>/abWarmStart.mat')
%
%   The study's runs deliberately did not feed the archive while it ran (all
%   arms used one frozen snapshot). Their final best rigs are valid CF3
%   solutions on the same normalisation table, and include the best rigs
%   found so far, so later warm starts should be able to use them. Each run's
%   BestChromosome is added with Source 'abstudy'; duplicates are skipped.
    S = load(resultFile, 'R');
    R = S.R;
    runs = R.runs(~isnan([R.runs.BestCost]));
    [~, stamp] = fileparts(fileparts(resultFile));
    n = 0;
    for nc = unique([runs.NumCams])
        rc = runs([runs.NumCams] == nc);
        key = conditionKey(nc, R.opts.CostFunc, R.opts.TargetType, R.opts.GridMode, R.opts.Spacing);
        for r = rc
            tag = sprintf('ABWarmStart/%s %s_%d rep%d', stamp, r.Method, r.NumSeeds, r.Repeat);
            n = n + updateSeedArchive(R.opts.Preset, key, r.BestChromosome, r.BestCost, ...
                'abstudy', tag, true);
        end
        pool = loadSeedArchive(R.opts.Preset, key);
        fprintf('%s: archive now %d solutions, best %.5f\n', key, numel(pool), pool(1).Cost);
    end
end
