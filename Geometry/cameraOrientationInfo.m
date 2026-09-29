function info = cameraOrientationInfo(cameraChromosome, numCams, degenerateTolDeg)
%CAMERAORIENTATIONINFO  Camera-frame axes, horizon tilt and inversion flags.
%
%   info = cameraOrientationInfo(chrom, numCams) decodes every camera's
%   orientation genes and reports, per camera, where the optical axis
%   points and how far the image is rolled away from level.
%
%   info = cameraOrientationInfo(chrom, numCams, degenerateTolDeg) sets the
%   tolerance (degrees, default 2) within which an optical axis counts as
%   vertical, i.e. the camera looks straight up or straight down and roll
%   relative to the horizon is undefined.
%
%   CONVENTIONS
%   The chromosome stores [x y z alpha beta gamma] per camera and
%   setupCameras builds R = eul2rotm([alpha beta gamma], "XYZ"), which is
%   the intrinsic sequence R = Rx(alpha)*Ry(beta)*Rz(gamma). The camera
%   frame follows Peter Corke's CentralCamera convention, so the columns of
%   R are the world directions of
%
%       R(:,1) = image RIGHT (+u)      R(:,2) = image DOWN (+v)
%       R(:,3) = OPTICAL AXIS (+z)
%
%   Because Rz is applied last, and Rz*[0;0;1] = [0;0;1], the optical axis
%   depends on alpha and beta ONLY. gamma is a pure roll about the optical
%   axis: changing it spins the image without moving where the camera
%   looks. That property is what uprightCameras exploits.
%
%   OUTPUT fields (each numCams rows)
%     Euler          numCams x 3  [alpha beta gamma] in radians
%     Right          numCams x 3  world direction of image +u
%     Down           numCams x 3  world direction of image +v
%     OpticalAxis    numCams x 3  world direction the camera looks along
%     Up             numCams x 3  world direction of image "up" (= -Down)
%     TiltDeg        numCams x 1  signed roll about the optical axis measured
%                                 from level, in (-180, 180]. 0 = horizon
%                                 level and upright, +/-180 = upside down.
%     Inverted       numCams x 1  logical, true when Up points below the
%                                 horizontal (|TiltDeg| > 90), i.e. the
%                                 camera perceives the world upside down.
%     ElevationDeg   numCams x 1  elevation of the optical axis above
%                                 horizontal; negative = looking downwards.
%     Degenerate     numCams x 1  logical, optical axis within
%                                 degenerateTolDeg of vertical, so TiltDeg
%                                 and Inverted are not meaningful.
%
%   See also: uprightCameras, snapChromosome, setupCameras.

    if nargin < 3 || isempty(degenerateTolDeg), degenerateTolDeg = 2; end

    worldUp = [0; 0; 1];
    sinTol  = sind(degenerateTolDeg);

    info.Euler        = zeros(numCams, 3);
    info.Right        = zeros(numCams, 3);
    info.Down         = zeros(numCams, 3);
    info.OpticalAxis  = zeros(numCams, 3);
    info.Up           = zeros(numCams, 3);
    info.TiltDeg      = zeros(numCams, 1);
    info.Inverted     = false(numCams, 1);
    info.ElevationDeg = zeros(numCams, 1);
    info.Degenerate   = false(numCams, 1);

    for c = 1:numCams
        idx = (c-1)*6 + 1;
        eul = cameraChromosome(idx+3:idx+5);
        R   = eul2rotm(eul(:).', "XYZ");

        xc = R(:,1);            % image right
        yc = R(:,2);            % image down
        zc = R(:,3);            % optical axis
        up = -yc;

        info.Euler(c,:)       = eul(:).';
        info.Right(c,:)       = xc.';
        info.Down(c,:)        = yc.';
        info.OpticalAxis(c,:) = zc.';
        info.Up(c,:)          = up.';
        info.ElevationDeg(c)  = asind(max(-1, min(1, zc(3))));

        % Level reference: the image-right direction the camera would have
        % if it were rolled so the horizon is horizontal in frame. It is
        % the horizontal direction perpendicular to the optical axis.
        h = cross(zc, worldUp);
        if norm(h) < sinTol
            % Optical axis (near) vertical: every roll leaves "up" in the
            % horizontal plane, so upright/inverted has no meaning here.
            info.Degenerate(c) = true;
            info.TiltDeg(c)    = NaN;
            info.Inverted(c)   = up(3) < 0;   % best effort, not meaningful
            continue;
        end

        rLevel = h / norm(h);
        % Signed roll from rLevel to the actual image-right, measured about
        % the optical axis with the right-hand rule.
        info.TiltDeg(c)  = atan2d(dot(cross(rLevel, xc), zc), dot(rLevel, xc));
        info.Inverted(c) = up(3) < 0;
    end
end
