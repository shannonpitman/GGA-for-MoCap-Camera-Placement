function d = runResultsDir(preset, numCams)
%RUNRESULTSDIR  Folder for per-run result files.
%   The lab preset keeps the historical Results/<N>Cams/; any other preset
%   gets Results/<preset>/<N>Cams/, so runs on different hardware never mix.
%   d = runResultsDir(preset) returns the preset's root (no <N>Cams).
    root = addProjectPaths();
    if nargin < 1 || isempty(preset) || strcmpi(preset, 'optitrack_lab')
        base = fullfile(root, 'Results');
    else
        base = fullfile(root, 'Results', preset);
    end
    if nargin < 2 || isempty(numCams)
        d = base;
    else
        d = fullfile(base, sprintf('%dCams', numCams));
    end
end
