function [chrom, report] = snapChromosome(chrom, varargin)
%SNAPCHROMOSOME  Snap a configuration to a mountable manufacturing grid.
%
%   [chrom, report] = snapChromosome(chrom, 'OrientationStepDeg', 15) rounds
%   every Euler gene to the nearest multiple of 15 degrees and reports how
%   far each camera actually moved.
%
%   WHY
%   The GA returns orientations to machine precision — 1.5708 rad, 0.7846
%   rad — and nobody mounts a camera to four decimal places. A bracket with
%   detents, a protractor, or a printed mounting plate realises discrete
%   angles. Snapping to the coarsest grid the cost can tolerate is what
%   makes an optimised layout physically reproducible, and the cost paid
%   for that is a result worth reporting.
%
%   Name-Value parameters
%     'OrientationStepDeg' Angular grid in degrees. Default 15. Use [] or 0
%                          to leave orientations alone.
%     'PositionStep'       Position grid in metres. Default [] (no snap).
%     'Genes'              Which orientation genes to snap, as a subset of
%                          [1 2 3] = [alpha beta gamma]. Default [1 2 3].
%                          Splitting them matters: alpha and beta steer the
%                          optical axis and so drive coverage, while gamma
%                          is a pure roll about that axis and is very nearly
%                          free. Pass [1 2] to price the pointing grid alone,
%                          or 3 to confirm the roll grid costs nothing.
%     'NumCams'            Camera count. Default numel(chrom)/6.
%     'VarMin', 'VarMax'   Optional 1 x 6*numCams search bounds. When given,
%                          the snapped chromosome is clamped to them, so a
%                          gene sitting exactly on a bound can never be
%                          rounded outside the feasible set.
%     'KeepUpright'        Re-run uprightCameras('Mode','flip') after
%                          snapping. Default true. This is a no-op whenever
%                          the step divides 180 degrees (1, 2, 5, 10, 15,
%                          30, 45, 90 all do), because a 180-degree flip
%                          then lands back on the grid; it matters only for
%                          steps such as 7 or 20 that do not.
%     'Verbose'            Print a per-camera table. Default false.
%
%   REPORT fields (per camera unless noted)
%     OrientationStepDeg, PositionStep      the grid actually applied
%     DeltaEulerDeg    numCams x 3  signed change in each Euler gene
%     GeodesicDeg      numCams x 1  true angle of the rotation taking the
%                                   original camera frame to the snapped
%                                   one — the honest "how much did this
%                                   camera turn" number, since Euler gene
%                                   deltas are not additive
%     AxisShiftDeg     numCams x 1  angle the OPTICAL AXIS moved. This is
%                                   the part that changes coverage; the
%                                   remainder of GeodesicDeg is roll.
%     RollShiftDeg     numCams x 1  change in roll about the optical axis,
%                                   wrapped into (-90, 90]. It is wrapped
%                                   because a 180-degree roll is free: it
%                                   maps the sensor rectangle onto itself.
%                                   Without the wrap, a camera that snapping
%                                   nudges across the inversion boundary and
%                                   KeepUpright then flips back reports a
%                                   ~170-degree rotation for what is
%                                   physically a small change of mount roll.
%                                   NaN for a camera whose optical axis is
%                                   vertical, where roll has no reference.
%     DeltaPos         numCams x 3  signed position change (m), zeros when
%                                   no position snap was requested
%     PosShift         numCams x 1  Euclidean position change (m)
%     GridResidualDeg  numCams x 3  distance from each ORIGINAL gene to the
%                                   nearest grid multiple, before snapping.
%                                   Small values mean the GA had already
%                                   landed near a mountable angle.
%     Clamped          logical, true if clamping to VarMin/VarMax bit
%     UprightFixed     number of cameras re-flipped by KeepUpright
%
%   See also: uprightCameras, cameraOrientationInfo, analyseConfiguration.

    p = inputParser;
    addParameter(p, 'OrientationStepDeg', 15,    @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'PositionStep',       [],    @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'Genes',              [1 2 3], @(x) isnumeric(x) && all(ismember(x, 1:3)));
    addParameter(p, 'NumCams',            [],    @isnumeric);
    addParameter(p, 'VarMin',             [],    @isnumeric);
    addParameter(p, 'VarMax',             [],    @isnumeric);
    addParameter(p, 'KeepUpright',        true,  @islogical);
    addParameter(p, 'Verbose',            false, @islogical);
    parse(p, varargin{:});
    opts = p.Results;

    if isempty(opts.NumCams)
        assert(mod(numel(chrom), 6) == 0, 'snapChromosome:BadLength', ...
            'Chromosome length %d is not a multiple of 6.', numel(chrom));
        numCams = numel(chrom) / 6;
    else
        numCams = opts.NumCams;
    end

    oriStepDeg = opts.OrientationStepDeg;
    if isempty(oriStepDeg), oriStepDeg = 0; end
    posStep = opts.PositionStep;
    if isempty(posStep), posStep = 0; end
    genes = unique(opts.Genes(:).');

    original = chrom;

    deltaEulerDeg   = zeros(numCams, 3);
    gridResidualDeg = zeros(numCams, 3);
    deltaPos        = zeros(numCams, 3);

    oriStepRad = deg2rad(oriStepDeg);

    for c = 1:numCams
        idx = (c-1)*6 + 1;

        if posStep > 0
            pos = chrom(idx:idx+2);
            chrom(idx:idx+2) = round(pos / posStep) * posStep;
        end

        if oriStepRad > 0 && ~isempty(genes)
            eul = chrom(idx+3:idx+5);
            % Distance from each original gene to its nearest grid multiple.
            gridResidualDeg(c,:) = rad2deg(abs(eul - round(eul / oriStepRad) * oriStepRad));
            eul(genes) = round(eul(genes) / oriStepRad) * oriStepRad;
            chrom(idx+3:idx+5) = eul;
        end
    end

    % Clamp before measuring displacement, so the reported deltas describe
    % the chromosome that is actually returned.
    clamped = false;
    if ~isempty(opts.VarMin) && ~isempty(opts.VarMax)
        clampedChrom = min(max(chrom, opts.VarMin), opts.VarMax);
        clamped = any(clampedChrom ~= chrom);
        chrom = clampedChrom;
    end

    upright = struct('NumChanged', 0);
    if opts.KeepUpright && oriStepRad > 0 && ~isempty(genes)
        [chrom, upright] = uprightCameras(chrom, 'NumCams', numCams, 'Mode', 'flip');
    end

    % Displacement measures, computed from the frames themselves rather
    % than from Euler differences (Euler gene deltas do not compose).
    geodesicDeg  = zeros(numCams, 1);
    axisShiftDeg = zeros(numCams, 1);
    posShift     = zeros(numCams, 1);

    tiltBefore = cameraOrientationInfo(original, numCams).TiltDeg;
    tiltAfter  = cameraOrientationInfo(chrom,    numCams).TiltDeg;
    % Wrap into (-90, 90]: a 180-degree roll leaves the field of view
    % unchanged, so it should not count as displacement.
    rollShiftDeg = mod(tiltAfter - tiltBefore + 90, 180) - 90;
    for c = 1:numCams
        idx = (c-1)*6 + 1;

        deltaPos(c,:) = chrom(idx:idx+2) - original(idx:idx+2);
        posShift(c)   = norm(deltaPos(c,:));

        deltaEulerDeg(c,:) = rad2deg(wrapAnglePi(chrom(idx+3:idx+5) - original(idx+3:idx+5)));

        R0 = eul2rotm(original(idx+3:idx+5), "XYZ");
        R1 = eul2rotm(chrom(idx+3:idx+5),    "XYZ");

        cosTheta = (trace(R0.' * R1) - 1) / 2;
        geodesicDeg(c) = acosd(max(-1, min(1, cosTheta)));

        axisShiftDeg(c) = acosd(max(-1, min(1, dot(R0(:,3), R1(:,3)))));
    end

    report.OrientationStepDeg = oriStepDeg;
    report.PositionStep       = posStep;
    report.Genes              = genes;
    report.DeltaEulerDeg      = deltaEulerDeg;
    report.GeodesicDeg        = geodesicDeg;
    report.AxisShiftDeg       = axisShiftDeg;
    report.RollShiftDeg       = rollShiftDeg;
    report.DeltaPos           = deltaPos;
    report.PosShift           = posShift;
    report.GridResidualDeg    = gridResidualDeg;
    report.Clamped            = clamped;
    report.UprightFixed       = upright.NumChanged;

    if opts.Verbose
        geneNames = {'alpha','beta','gamma'};
        fprintf('\n  snapChromosome — orientation %g deg on [%s]', ...
            oriStepDeg, strjoin(geneNames(genes), ' '));
        if posStep > 0, fprintf(', position %g m', posStep); end
        fprintf('\n');
        fprintf('  %-5s  %-9s  %-10s  %-10s  %-10s  %s\n', ...
            'Cam', 'Rotated', 'AxisMoved', 'RollMoved', 'Moved(m)', 'Grid residual a/b/g (deg)');
        fprintf('  %s\n', repmat('-', 1, 86));
        for c = 1:numCams
            fprintf('  %-5d  %8.2f   %9.2f   %9.2f   %9.3f   %5.2f %5.2f %5.2f\n', ...
                c, geodesicDeg(c), axisShiftDeg(c), rollShiftDeg(c), posShift(c), ...
                gridResidualDeg(c,1), gridResidualDeg(c,2), gridResidualDeg(c,3));
        end
        fprintf(['  max optical-axis shift %.2f deg (this is what moves coverage),\n' ...
                 '  max roll shift %.2f deg, max total rotation as applied %.2f deg\n\n'], ...
            max(axisShiftDeg), max(abs(rollShiftDeg), [], 'omitnan'), max(geodesicDeg));
    end
end
