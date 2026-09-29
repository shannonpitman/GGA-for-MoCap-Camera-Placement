function r = rotvecWrap(r)
%ROTVECWRAP  Unique rotation vector with |r| <= pi (same rotation).
%   A rotation by theta about k equals a rotation by theta - 2*pi about k,
%   so any r with |r| > pi is mapped back into the pi-ball. Replaces angle
%   clamping for the exponential-map parameterisation.
    theta = norm(r);
    if theta <= pi
        return;
    end
    k = r / theta;
    theta = mod(theta + pi, 2*pi) - pi;  % into (-pi, pi]
    r = k * theta;
end
