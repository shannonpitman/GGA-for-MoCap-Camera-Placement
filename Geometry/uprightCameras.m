function [chrom, report] = uprightCameras(chrom, varargin)
%UPRIGHTCAMERAS  Roll cameras about their optical axes so none is inverted.
%
%   [chrom, report] = uprightCameras(chrom) returns a chromosome in which
%   no camera perceives the world upside down. Only each camera's roll
%   about its own optical axis changes.
%
%   WHY THIS EXISTS
%   The GA never constrains roll. Its cost function only asks which target
%   points fall inside each image, and that is (all but exactly) invariant
%   to spinning a camera about its own optical axis, so the search is free
%   to return perfectly good configurations in which cameras hang upside
%   down. That is fine numerically and useless physically: an inverted
%   camera has to be mounted inverted, its status LEDs and cable gland face
%   the wrong way, and Motive's 2-D views are unreadable during calibration.
%   This function removes that artefact after the fact rather than
%   constraining the GA, so previously logged runs stay valid.
%
%   WHY IT IS (ESSENTIALLY) FREE
%   The roll is applied as R -> R*Rz(delta), a rotation about the camera's
%   own +z — the optical axis — so it rolls the image without moving the
%   optical axis at all (identical to the old gamma -> gamma + delta on XYZ
%   Euler genes). In 'flip' mode delta = pi, a 180-degree image rotation. A 180-degree rotation
%   maps the sensor rectangle onto itself, so the set of directions inside
%   the field of view is unchanged and coverage is preserved to floating
%   point (the only exception is a target that lands exactly on the far
%   pixel edge, since the visibility test is u,v in [1, W] about a
%   principal point at W/2 — a measure-zero case worth nothing).
%
%   MODES
%     'flip'  (default) Roll inverted cameras by 180 degrees only.
%             Coverage-neutral. Leaves the residual roll away from level
%             untouched, so cameras stay tilted exactly as the GA left them,
%             just no longer inverted.
%     'level' Roll every camera so its roll about the optical axis is
%             zero, i.e. the horizon is level in frame. This is NOT
%             coverage-neutral: the sensor is 1280 x 1024, so rolling by
%             anything other than a multiple of 180 degrees re-orients a
%             non-square field of view and the cost changes slightly. Use
%             it to quantify that penalty, not as the default fix.
%
%   Name-Value parameters
%     'Mode'             'flip' (default) or 'level'
%     'DegenerateTolDeg' Optical axes within this angle of vertical are
%                        left untouched, because "upside down" is undefined
%                        for a camera looking straight up or down.
%                        Default 2 degrees.
%     'NumCams'          Camera count. Default numel(chrom)/6.
%     'Verbose'          Print a per-camera table. Default false.
%
%   REPORT fields
%     NumCams, Mode, Changed (logical per camera), NumChanged,
%     DeltaGammaDeg, TiltBeforeDeg, TiltAfterDeg, InvertedBefore,
%     InvertedAfter, Degenerate, ElevationDeg
%
%   See also: cameraOrientationInfo, snapChromosome, analyseConfiguration.

    p = inputParser;
    addParameter(p, 'Mode',             'flip', @(s) any(strcmpi(s, {'flip','level'})));
    addParameter(p, 'DegenerateTolDeg', 2,      @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'NumCams',          [],     @isnumeric);
    addParameter(p, 'Verbose',          false,  @islogical);
    parse(p, varargin{:});
    opts = p.Results;

    if isempty(opts.NumCams)
        assert(mod(numel(chrom), 6) == 0, 'uprightCameras:BadLength', ...
            'Chromosome length %d is not a multiple of 6.', numel(chrom));
        numCams = numel(chrom) / 6;
    else
        numCams = opts.NumCams;
    end

    before = cameraOrientationInfo(chrom, numCams, opts.DegenerateTolDeg);

    deltaDeg = zeros(numCams, 1);
    for c = 1:numCams
        if before.Degenerate(c)
            continue;   % roll relative to the horizon is undefined here
        end
        switch lower(opts.Mode)
            case 'flip'
                if before.Inverted(c)
                    deltaDeg(c) = 180;
                end
            case 'level'
                deltaDeg(c) = -before.TiltDeg(c);
        end
    end

    for c = 1:numCams
        if deltaDeg(c) == 0, continue; end
        rIdx = (c-1)*6 + (4:6);
        d = deg2rad(deltaDeg(c));
        Rroll = [cos(d) -sin(d) 0; sin(d) cos(d) 0; 0 0 1];   % about camera +z
        chrom(rIdx) = rotmToGenes(genesToRotm(chrom(rIdx)) * Rroll);
    end

    after = cameraOrientationInfo(chrom, numCams, opts.DegenerateTolDeg);

    report.NumCams        = numCams;
    report.Mode           = lower(opts.Mode);
    report.Changed        = deltaDeg ~= 0;
    report.NumChanged     = nnz(deltaDeg ~= 0);
    report.DeltaGammaDeg  = deltaDeg;
    report.TiltBeforeDeg  = before.TiltDeg;
    report.TiltAfterDeg   = after.TiltDeg;
    report.InvertedBefore = before.Inverted;
    report.InvertedAfter  = after.Inverted;
    report.Degenerate     = before.Degenerate;
    report.ElevationDeg   = before.ElevationDeg;

    if opts.Verbose
        fprintf('\n  uprightCameras (mode: %s) — %d of %d cameras adjusted\n', ...
            report.Mode, report.NumChanged, numCams);
        fprintf('  %-5s  %-10s  %-10s  %-10s  %-9s  %s\n', ...
            'Cam', 'Tilt(before)', 'Tilt(after)', 'dGamma', 'Elev', 'Note');
        fprintf('  %s\n', repmat('-', 1, 68));
        for c = 1:numCams
            if report.Degenerate(c)
                note = 'vertical axis — left alone';
            elseif report.InvertedBefore(c)
                note = 'was INVERTED';
            else
                note = '';
            end
            fprintf('  %-5d  %10.1f  %10.1f  %10.1f  %9.1f  %s\n', ...
                c, report.TiltBeforeDeg(c), report.TiltAfterDeg(c), ...
                report.DeltaGammaDeg(c), report.ElevationDeg(c), note);
        end
        fprintf('\n');
    end
end
