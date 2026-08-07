function projectRoot = addProjectPaths()
% addProjectPaths  Add every code subfolder of this project to the MATLAB path.
%
%   Call this once at the start of a session (or from any entry-point script
%   such as runCameraOptimiser, batchRunGA, plotGARuns, analyseConfiguration)
%   so that functions in Run/, ParameterTesting/, GA_Core/, CostFunctions/,
%   Geometry/, Setup/, Plotting/, Analysis/ and Sensitivity/ resolve
%   regardless of which folder MATLAB is currently in.
%
%   projectRoot = addProjectPaths() also RETURNS the absolute path of the
%   project root. Use that return value whenever you need to build a path to
%   Results/ or figures/ — do NOT use fileparts(mfilename('fullpath')) from
%   inside a subfolder, because that yields the subfolder, not the root.
%   Every script in this project uses the return value, so files can be moved
%   between subfolders without breaking their Results/ paths.
%
%   This file must stay at the project root: it is what defines where the
%   root is.
%
%   Excludes _Archive/ on purpose — anything in there is intentionally
%   off the active path.
%
%   You normally do not call this yourself: START_HERE (in the workspace
%   root) calls it for you, after putting the RVC3 toolbox on the path.

    projectRoot = fileparts(mfilename('fullpath'));

    codeSubfolders = {'Run',      'ParameterTesting', ...
                      'GA_Core',  'CostFunctions',    'Geometry', ...
                      'Setup',    'Plotting',         'Analysis', ...
                      'Sensitivity'};

    for i = 1:numel(codeSubfolders)
        p = fullfile(projectRoot, codeSubfolders{i});
        if isfolder(p)
            addpath(p);
        end
    end

    % Also put the project root on the path so addProjectPaths itself stays
    % reachable after a cd.
    addpath(projectRoot);
end
