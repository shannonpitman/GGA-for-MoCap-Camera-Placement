function f = normTableFile(preset)
%NORMTABLEFILE  Utopia/nadir table for a runConfig preset.
%   The lab preset keeps the historical Results/normTable.mat; every other
%   preset gets its own file, because the constants depend on hardware,
%   workspace and mount model.
    root = addProjectPaths();
    if nargin < 1 || isempty(preset) || strcmpi(preset, 'optitrack_lab')
        f = fullfile(root, 'Results', 'normTable.mat');
    else
        f = fullfile(root, 'Results', sprintf('normTable_%s.mat', preset));
    end
end
