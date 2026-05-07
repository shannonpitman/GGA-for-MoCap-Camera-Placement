function projectRoot = addProjectPaths()
% addProjectPaths  Add every code subfolder of this project to the MATLAB path.
%
%   Call this once at the start of a session (or from any entry-point script
%   such as runCameraOptimiser, batchRunGA, plotGARuns, analyseConfiguration)
%   so that functions in GA_Core/, CostFunctions/, Geometry/, Setup/,
%   Plotting/ and Analysis/ resolve regardless of which folder MATLAB is in.
%
%   Returns the absolute path of the project root, useful for building
%   paths to the Results/ and figures/ folders.
%
%   Excludes Unused/ on purpose — anything in there is intentionally
%   off the active path.

    projectRoot = fileparts(mfilename('fullpath'));

    codeSubfolders = {'GA_Core', 'CostFunctions', 'Geometry', ...
                      'Setup',   'Plotting',      'Analysis', ...
                      'Sensitivity'};

    for i = 1:numel(codeSubfolders)
        p = fullfile(projectRoot, codeSubfolders{i});
        if isfolder(p)
            addpath(p);
        end
    end

    % Also put the project root on the path so the entry-point scripts
    % (runCameraOptimiser, batchRunGA, plotGARuns, analyseConfiguration,
    % ConfigAnalyser, viewGALog, app1-Shannons-PC) are reachable.
    addpath(projectRoot);
end
