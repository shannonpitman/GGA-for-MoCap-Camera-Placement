function tests = testWarmSeeds
% Seed selection must inject exactly the requested number of rows, respect
% the diversity threshold for 'top_n', and reproduce the old scheme for
% 'best_mutant'. The archive must deduplicate.
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    addProjectPaths();
    rng(1);
    nc = 3;
    problem.nVar = 6*nc;
    problem.VarMin = repmat([-5 -4.5 0 -pi -pi -pi], 1, nc);
    problem.VarMax = repmat([ 5  4.5 4.8 pi pi pi], 1, nc);
    base = [ -3 -3 2 0 0 0,  3 -3 2 0 0 0,  0 3 2 0 0 0 ];
    pool = struct('Chromosome', {}, 'Cost', {});
    % 6 distinct layouts (shifted 1 m each), each with a near-copy 5 cm away
    for k = 1:6
        c = base;  c(1:6:end) = c(1:6:end) + (k-1);
        pool(end+1) = struct('Chromosome', c, 'Cost', k); %#ok<AGROW>
        c2 = c;  c2(2:6:end) = c2(2:6:end) + 0.05;
        pool(end+1) = struct('Chromosome', c2, 'Cost', k + 0.5); %#ok<AGROW>
    end
    tc.TestData.problem = problem;
    tc.TestData.pool = pool(randperm(numel(pool)));
end

function testTopNDiverseFirst(tc)
    [S, info] = selectWarmSeeds(tc.TestData.pool, 6, 'top_n', 0.25, tc.TestData.problem);
    verifySize(tc, S, [6 18]);
    verifyEqual(tc, info.nDistinct, 6);
    verifyEqual(tc, info.Costs.', 1:6);          % the near-copies were skipped
end

function testTopNFillsWithRejected(tc)
    [S, info] = selectWarmSeeds(tc.TestData.pool, 9, 'top_n', 0.25, tc.TestData.problem);
    verifySize(tc, S, [9 18]);
    verifyEqual(tc, [info.nDistinct info.nFilled info.nMutant], [6 3 0]);
    verifyEqual(tc, info.Costs(7:9).', [1.5 2.5 3.5]);
end

function testTopNNoDistanceIsPureRanking(tc)
    [~, info] = selectWarmSeeds(tc.TestData.pool, 4, 'top_n', 0, tc.TestData.problem);
    verifyEqual(tc, info.Costs.', [1 1.5 2 2.5]);
end

function testTopNPoolTooSmallPadsWithMutants(tc)
    [S, info] = selectWarmSeeds(tc.TestData.pool, 15, 'top_n', 0.25, tc.TestData.problem);
    verifySize(tc, S, [15 18]);
    verifyEqual(tc, info.nMutant, 3);
    M = S(13:15, :);   % the mutants (the fixture's own layouts run past x = 5)
    verifyTrue(tc, all(M >= tc.TestData.problem.VarMin & M <= tc.TestData.problem.VarMax, 'all'));
end

function testBestMutant(tc)
    [S, info] = selectWarmSeeds(tc.TestData.pool, 10, 'best_mutant', 0.25, tc.TestData.problem);
    verifySize(tc, S, [10 18]);
    best = tc.TestData.pool([tc.TestData.pool.Cost] == 1).Chromosome;
    verifyEqual(tc, S(1, :), best);
    verifyEqual(tc, info.nMutant, 9);
    verifyGreaterThan(tc, min(vecnorm(S(2:end, :) - best, 2, 2)), 0);
end

function testLayoutDistanceIgnoresSlotOrder(tc)
    a = [ -3 -3 2 0 0 0,  3 -3 2 0.1 0 0,  0 3 2 0 0 0 ];
    b = [a(13:18), a(1:6), a(7:12)];
    verifyEqual(tc, layoutDistance(a, b), 0, 'AbsTol', 1e-12);
    c = a;  c(1) = c(1) + 0.3;
    verifyEqual(tc, layoutDistance(a, c), 0.1, 'AbsTol', 1e-12);
end

function testArchiveDeduplicates(tc)
    preset = sprintf('test_%s', char(java.util.UUID.randomUUID));
    key = conditionKey(3, 3, 1, 1, 1.0);
    cleanup = onCleanup(@() rmdir(fileparts(seedArchiveFile(preset, key)), 's')); %#ok<NASGU>
    C = vertcat(tc.TestData.pool(1:3).Chromosome);  J = [tc.TestData.pool(1:3).Cost].';
    verifyEqual(tc, updateSeedArchive(preset, key, C, J, 'best', 'a.mat', false), 3);
    verifyEqual(tc, updateSeedArchive(preset, key, C, J, 'best', 'a.mat', false), 0);
    pool = loadSeedArchive(preset, key);
    verifyEqual(tc, numel(pool), 3);
    verifyTrue(tc, issorted([pool.Cost]));
end
