function r = aimRotvec(direction)
%AIMROTVEC  Rotation vector that turns the camera optical axis (+z) onto
%   a world direction: axis = z x d, angle = acos(z . d). This is the
%   axis-angle used to aim cameras at section centres; it is stored
%   directly as the orientation genes (no quaternion/Euler conversion).
    d = direction(:).' / norm(direction);
    z = [0 0 1];
    axis = cross(z, d);
    s = norm(axis);
    angle = atan2(s, dot(z, d));
    if s < 1e-12
        if dot(z, d) > 0
            r = [0 0 0];                 % already looking along +z
        else
            r = [pi 0 0];                % looking straight down: 180 deg about x
        end
        return;
    end
    r = axis / s * angle;
end
