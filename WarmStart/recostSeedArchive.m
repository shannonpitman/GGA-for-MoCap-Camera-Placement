function recostSeedArchive(preset, varargin)
%RECOSTSEEDARCHIVE  Re-score archived CF3 rigs on the current normalisation table.
%
%   recostSeedArchive('optitrack_lab')
%
%   CF3 costs depend on the utopia/nadir table, so after the table changes
%   (buildNormFromArchive) the stored CF3 costs in the seed archive are on
%   the old scale and would rank warm-start seeds wrongly. This re-evaluates
%   every CF3 archive entry with the current table and re-sorts. CF1/CF2
%   costs are raw and unaffected. The old files are copied to
%   Results/SeedArchive/<preset>/_backup_<stamp>/ first.
    root = addProjectPaths(); %#ok<NASGU>
    cfg = runConfig(preset);
    archDir = fileparts(seedArchiveFile(preset, 'x'));
    files = dir(fullfile(archDir, '*C_CF3_TT*_GM*_sp*.mat'));
    if isempty(files)
        fprintf('No CF3 archives for %s.\n', preset);
        return;
    end
    bak = fullfile(archDir, ['_backup_' char(string(datetime('now'), 'yyyyMMdd_HHmmss'))]);
    mkdir(bak);
    for i = 1:numel(files)
        copyfile(fullfile(archDir, files(i).name), bak);
        k = sscanf(files(i).name, '%dC_CF%d_TT%d_GM%d_sp%f.mat');
        S = load(fullfile(archDir, files(i).name), 'archive');
        archive = S.archive;
        specs = buildRunSpecs(cfg, k(1), 3, k(3), k(4), k(5));
        C = vertcat(archive.Chromosome);
        J = zeros(size(C, 1), 1);
        parfor j = 1:size(C, 1)
            J(j) = combinedCostFunction(C(j, :), specs);
        end
        oldBest = archive(1).Cost;
        for j = 1:numel(archive), archive(j).Cost = J(j); end
        [~, order] = sort([archive.Cost]);
        archive = archive(order); %#ok<NASGU>
        save(fullfile(archDir, files(i).name), 'archive');
        fprintf('  %s: %d rigs, best %.5f (was %.5f)\n', files(i).name, numel(J), min(J), oldBest);
    end
    fprintf('Old archives copied to %s\n', bak);
end
