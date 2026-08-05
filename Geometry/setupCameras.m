function [cameras, cameraCenters] = setupCameras(cameraChromosome, numCams, resolution, focalLength, focalLengthWide, PrincipalPoint, pixelSize)
%SETUPCAMERAS  Initialise CentralCamera objects from a chromosome.
%
% [cameras, cameraCenters] = setupCameras(chrom, n, res, f, fw, c0, rho)
%
% pixelSize (rho) is OPTIONAL but strongly recommended. If omitted the
% function emits a one-shot warning and falls back to Peter Corke's
% CentralCamera default (10 um), which gives an unphysically wide
% simulated FOV (~99° narrow, ~122° wide) for the OptiTrack focals
% (5.5 mm / 3.5 mm) and resolution (1280 x 1024). With pixelSize =
% specs.PixelSize (4.57 µm derived from datasheet) the simulated FOV
% matches physical (~56° narrow, ~82° wide).
    cameras = cell(numCams, 1);
    cameraCenters = zeros(3, numCams);

    % Optional pixel size — fall back to Corke default with a warning so
    % legacy callers do not silently use the wrong FOV.
    if nargin < 7 || isempty(pixelSize)
        persistent warned
        if isempty(warned)
            warning('setupCameras:NoPixelSize', ...
                ['pixelSize (rho) not supplied; CentralCamera will use ' ...
                 'its default 10 um. Simulated FOV will be wider than ' ...
                 'physical. Pass specs.PixelSize through to avoid this.']);
            warned = true;
        end
        pixelSize = [];   % signal "do not pass to CentralCamera"
    end

    for i = 1:numCams
        chromStartIdx = (i-1)*6 + 1;
        chromEndIdx   = i*6;
        camPosition    = cameraChromosome(chromStartIdx:chromStartIdx+2);
        camOrientation = cameraChromosome(chromEndIdx-2:chromEndIdx);

        T = se3(eul2rotm(camOrientation, "XYZ"), camPosition);
        if mod(i,2) == 0
            f = focalLengthWide;
        else
            f = focalLength;
        end

        if isempty(pixelSize)
            cameras{i} = CentralCamera( ...
                name=sprintf('cam%d', i), ...
                resolution=resolution, focal=f, pose=T, ...
                center=PrincipalPoint);
        else
            cameras{i} = CentralCamera( ...
                name=sprintf('cam%d', i), ...
                resolution=resolution, focal=f, pose=T, ...
                center=PrincipalPoint, pixel=pixelSize);
        end
        cameraCenters(:, i) = camPosition(:);
    end
end