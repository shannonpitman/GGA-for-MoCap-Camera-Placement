function thesisFigureSet(varargin)
%THESISFIGURESET  Regenerate the whole dissertation figure set as PNG + PDF.
%
%   thesisFigureSet()
%   thesisFigureSet('Only', 'config')      % just the GA-vs-ad-hoc figures
%   thesisFigureSet('CostFunctions', 3)    % CF3 only (faster)
%
%   Every figure is written twice — <name>.pdf for LaTeX and <name>.png at
%   300 dpi for Word and slides — into <projectRoot>/figures/.
%
%   Naming: <type>_<CF>_<TT>_<GM>_<spacing>_<cameras>, so a figure's
%   instance is readable off the filename. Figures that are inherently a
%   comparison across one of those factors drop that factor from the name.
%
%   SETS
%     'process'  convergence, population diversity, cost box plots,
%                computation time, warm-vs-cold — the GA behaviour figures
%     'config'   coverage heat-map, cost field, camera poses + FOV,
%                pairwise baseline angles, per-camera angle sensitivity,
%                cost composition — the GA-vs-ad-hoc figures at k = 7
%     'all'      both (default)
%
%   Name-Value parameters:
%     'Only'          'all' | 'process' | 'config'.   Default 'all'
%     'CameraRange'   camera counts.                  Default [6 7]
%     'CostFunctions' cost functions for the process set. Default [1 2 3]
%     'TargetTypes'   1 UAV, 2 UGV.                   Default [1 2]
%     'GridModes'     1 uniform, 2 normal.            Default [1 2]
%     'Spacing'       grid spacing [m].               Default 1.0
%     'RepCams'       camera count for per-run figures. Default 7
%     'OutDir'        Default <projectRoot>/figures
%
%   This supersedes plotGARuns.m for chapter output: plotGARuns sweeps
%   spacings that were never run and emits PDFs only. Keep plotGARuns for
%   a raw dump of everything the log happens to contain.

    projectRoot = addProjectPaths();

    p = inputParser;
    addParameter(p, 'Only',          'all',   @ischar);
    addParameter(p, 'CameraRange',   [6 7],   @isnumeric);
    addParameter(p, 'CostFunctions', [1 2 3], @isnumeric);
    addParameter(p, 'TargetTypes',   [1 2],   @isnumeric);
    addParameter(p, 'GridModes',     [1 2],   @isnumeric);
    addParameter(p, 'Spacing',       1.0,     @isnumeric);
    addParameter(p, 'RepCams',       7,       @isnumeric);
    addParameter(p, 'OutDir',        fullfile(projectRoot,'figures'), @ischar);
    parse(p, varargin{:});
    o = p.Results;

    if ~isfolder(o.OutDir), mkdir(o.OutDir); end

    logFile = fullfile(projectRoot, 'Results', 'Logs', 'GGA_RunsLog.mat');
    runDir  = fullfile(projectRoot, 'Results');
    common  = {'LogFile', logFile};

    cfTag = {'CF1', 'CF2', 'CF3'};
    ttTag = {'UAV', 'UGV'};
    gmTag = {'GM1uniform', 'GM2normal'};
    spTag = sprintf('sp%.0fcm', o.Spacing*100);

    doProcess = any(strcmpi(o.Only, {'all', 'process'}));
    doConfig  = any(strcmpi(o.Only, {'all', 'config'}));

    nOK = 0; nFail = 0;
    failures = {};

    % Every figure goes through this so one bad cell cannot abort the set.
    function attempt(label, fn)
        try
            fn();
            nOK = nOK + 1;
        catch ME
            nFail = nFail + 1;
            failures{end+1} = sprintf('%s  ->  %s', label, ME.message); %#ok<AGROW>
            fprintf(2, '  SKIPPED %s: %s\n', label, ME.message);
        end
        close all;
    end

    %% ==================================================================
    if doProcess
    %% ==================================================================

    section('CONVERGENCE');
    for cf = o.CostFunctions(:)'
        for tt = o.TargetTypes(:)'
            for gm = o.GridModes(:)'
                for nc = o.CameraRange(:)'
                    tag = sprintf('convergence_%s_%s_%s_%s_%dC', ...
                        cfTag{cf}, ttTag{tt}, gmTag{gm}, spTag, nc);
                    attempt(tag, @() plotGA_Convergence(common{:}, ...
                        'CostFunction', cf, 'TargetType', tt, 'GridMode', gm, ...
                        'Spacing', o.Spacing, 'NumCameras', nc, ...
                        'RunDir', runDir, 'ShowTopTen', true, ...
                        'ShowAvgCost', false, 'LogScale', true, ...
                        'SaveAs', fullfile(o.OutDir, tag)));
                end
            end
        end
    end

    section('POPULATION DIVERSITY');
    for cf = o.CostFunctions(:)'
        for tt = o.TargetTypes(:)'
            for gm = o.GridModes(:)'
                tag = sprintf('diversity_%s_%s_%s_%s_%dC', ...
                    cfTag{cf}, ttTag{tt}, gmTag{gm}, spTag, o.RepCams);
                attempt(tag, @() plotGA_PopulationDiversity(common{:}, ...
                    'CostFunction', cf, 'TargetType', tt, 'GridMode', gm, ...
                    'Spacing', o.Spacing, 'NumCameras', o.RepCams, ...
                    'RunDir', runDir, 'OverlayCost', true, ...
                    'SaveAs', fullfile(o.OutDir, tag)));
            end
        end
    end

    section('COST BOX PLOTS (cost vs camera count)');
    for cf = o.CostFunctions(:)'
        for gm = o.GridModes(:)'
            tag = sprintf('costbox_byTargetType_%s_%s', gmTag{gm}, spTag);
            attempt(tag, @() plotGA_CostBoxPlots(common{:}, ...
                'CostFunction', cf, 'GridMode', gm, 'Spacing', o.Spacing, ...
                'SplitBy', 'TargetType', ...
                'SaveAs', fullfile(o.OutDir, tag)));
        end
        for tt = o.TargetTypes(:)'
            tag = sprintf('costbox_byGridMode_%s_%s', ttTag{tt}, spTag);
            attempt(tag, @() plotGA_CostBoxPlots(common{:}, ...
                'CostFunction', cf, 'TargetType', tt, 'Spacing', o.Spacing, ...
                'SplitBy', 'GridMode', ...
                'SaveAs', fullfile(o.OutDir, tag)));
        end
    end

    section('COMPUTATION TIME');
    for gm = o.GridModes(:)'
        tag = sprintf('comptime_byTargetType_%s_%s', gmTag{gm}, spTag);
        attempt(tag, @() plotGA_ComputationTime(common{:}, ...
            'GridMode', gm, 'Spacing', o.Spacing, 'SplitBy', 'TargetType', ...
            'SaveAs', fullfile(o.OutDir, tag)));

        tag = sprintf('comptime_byCostFunction_%s_%s', gmTag{gm}, spTag);
        attempt(tag, @() plotGA_ComputationTime(common{:}, ...
            'GridMode', gm, 'Spacing', o.Spacing, 'SplitBy', 'CostFunction', ...
            'SaveAs', fullfile(o.OutDir, tag)));
    end

    section('WARM-START vs COLD-START');
    for cf = o.CostFunctions(:)'
        for tt = o.TargetTypes(:)'
            for gm = o.GridModes(:)'
                tag = sprintf('warmcold_%s_%s_%s_%s', ...
                    cfTag{cf}, ttTag{tt}, gmTag{gm}, spTag);
                % No ad-hoc overlay here. This figure is a within-GA
                % comparison, and after the utopia/nadir renormalisation
                % the ad-hoc rig scores 10-60x the GA — plotting it would
                % flatten both boxes to a line. The GA-vs-ad-hoc
                % comparison has its own figures.
                attempt(tag, @() plotGA_WarmColdEffect(common{:}, ...
                    'CostFunction', cf, 'TargetType', tt, 'GridMode', gm, ...
                    'Spacing', o.Spacing, 'OptiTrackOverlay', false, ...
                    'SaveAs', fullfile(o.OutDir, tag)));
            end
        end
    end

    end % doProcess

    %% ==================================================================
    if doConfig
    %% ==================================================================
    cfg = [common, {'CostFunction', 3, 'NumCameras', o.RepCams, ...
                    'Spacing', o.Spacing}];

    section('COST COMPOSITION (GA vs ad-hoc)');
    attempt('costcomponents', @() plotCostComponents_GAvsOptiTrack( ...
        'NumCameras', o.RepCams, 'Spacing', o.Spacing, 'LogFile', logFile, ...
        'SaveAs', fullfile(o.OutDir, sprintf('costcomponents_GAvsAdhoc_%s_%dC', ...
                                             spTag, o.RepCams))));

    section('GA vs AD-HOC, PER INSTANCE');
    for tt = o.TargetTypes(:)'
        for gm = o.GridModes(:)'
            inst = sprintf('%s_%s_%s_%dC', ttTag{tt}, gmTag{gm}, spTag, o.RepCams);
            here = [cfg, {'TargetType', tt, 'GridMode', gm}];

            tag = ['coverage_GAvsAdhoc_' inst];
            attempt(tag, @() plotHeatmap_GAvsOptiTrack(here{:}, ...
                'SaveAs', fullfile(o.OutDir, tag)));

            tag = ['costfield_GAvsAdhoc_' inst];
            attempt(tag, @() plotCostField_GAvsOptiTrack(here{:}, ...
                'SaveAs', fullfile(o.OutDir, tag)));

            tag = ['configfov_GAvsAdhoc_' inst];
            attempt(tag, @() plotConfigFOV_GAvsOptiTrack(here{:}, ...
                'SaveAs', fullfile(o.OutDir, tag)));

            tag = ['baselineangles_GAvsAdhoc_' inst];
            attempt(tag, @() plotBaselineAngles_GAvsOptiTrack(here{:}, ...
                'SaveAs', fullfile(o.OutDir, tag)));

            tag = ['anglesensitivity_cam1_' inst];
            attempt(tag, @() plotAngleSensitivity(here{:}, ...
                'CameraIndex', 1, 'SweepRange', 45, 'SweepStep', 3, ...
                'ShowAllCams', false, ...
                'SaveAs', fullfile(o.OutDir, tag)));
        end
    end

    section('GRID-MODE COVERAGE COMPARISON');
    tag = sprintf('coverage_uniformVsNormal_UAV_%s_%dC', spTag, o.RepCams);
    attempt(tag, @() plotCoverageHeatmap(common{:}, ...
        'NumCameras', o.RepCams, 'CostFunction', 3, ...
        'SaveAs', fullfile(o.OutDir, tag)));

    end % doConfig

    %% ==================================================================
    section('DONE');
    fprintf('  %d figure(s) written, %d skipped.\n', nOK, nFail);
    for i = 1:numel(failures)
        fprintf('    ! %s\n', failures{i});
    end
    pngs = dir(fullfile(o.OutDir, '*.png'));
    fprintf('  %d PNG files now in %s\n', numel(pngs), o.OutDir);
end


function section(txt)
    fprintf('\n===== %s =====\n', txt);
end
