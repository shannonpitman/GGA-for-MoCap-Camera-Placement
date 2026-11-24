function f = seedArchiveFile(preset, condKey)
%SEEDARCHIVEFILE  Results/SeedArchive/<preset>/<condKey>.mat
%   One archive per hardware/workspace preset and condition: costs are only
%   comparable within the same preset, cost function, target and grid.
    if isempty(preset), preset = 'optitrack_lab'; end
    root = addProjectPaths();
    d = fullfile(root, 'Results', 'SeedArchive', preset);
    old = fullfile(root, 'Results', 'Unconstrained', 'SeedArchive', preset);
    if ~isfolder(d) && isfolder(old), d = old; end   % superseded any-tilt OpenMV sweeps
    f = fullfile(d, [condKey '.mat']);
end
