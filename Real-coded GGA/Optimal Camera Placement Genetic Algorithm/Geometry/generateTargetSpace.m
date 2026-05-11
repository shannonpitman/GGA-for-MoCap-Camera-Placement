function targetSpace = generateTargetSpace(volume, targetMode, spacing)
% GENERATETARGETSPACE  Build the Nx3 grid of [x y z] target points.
%
% targetSpace = generateTargetSpace(volume, targetMode, spacing)
%
%   volume     : 3x2 matrix of [xLo xHi; yLo yHi; zLo zHi]   (metres)
%   targetMode : 1 = uniform grid
%                2 = normally-concentrated grid (inverse-CDF warping)
%   spacing    : scalar (isotropic, applied to x, y, z) or
%                1x3 / 3x1 vector [sx sy sz] for anisotropic spacing.
%
% The anisotropic form is used by the UGV target space so that x-y can be
% coarse (matching the UAV x-y resolution for thesis consistency) while z
% keeps enough layers to sample the floor-slab thickness — typically
% z = [0, 0.25, 0.5] m  ->  3 layers on a 0.5 m slab.

% Normalise spacing to a 1x3 vector
if isscalar(spacing)
    spacing = [spacing spacing spacing];
else
    spacing = spacing(:).';
    assert(numel(spacing) == 3, ...
        'generateTargetSpace:BadSpacing', ...
        'spacing must be a scalar or a 1x3 / 3x1 vector. Got %d elements.', ...
        numel(spacing));
end
assert(all(spacing > 0), 'generateTargetSpace:BadSpacing', ...
    'All spacing components must be > 0. Got [%g %g %g].', spacing);

% Generate the marker coordinates along each axis
x_marker = volume(1,1):spacing(1):volume(1,2);
y_marker = volume(2,1):spacing(2):volume(2,2);
z_marker = volume(3,1):spacing(3):volume(3,2);

% Guard against degenerate axes: if rounding drops the upper bound, force
% inclusion. Also ensure at least one z-layer exists when the slab is
% thinner than the spacing.
if isempty(x_marker), x_marker = volume(1,1); end
if isempty(y_marker), y_marker = volume(2,1); end
if isempty(z_marker), z_marker = volume(3,1); end

switch targetMode
    case 1
        % Uniform Grid
        [X, Y, Z] = meshgrid(x_marker, y_marker, z_marker);
        targetSpace = [X(:), Y(:), Z(:)];

    case 2
        % Normally-concentrated grid via inverse-CDF warping. Point counts
        % per axis are taken from the uniform grid to keep the discretisation
        % comparable.
        Nx = numel(x_marker);
        Ny = numel(y_marker);
        Nz = numel(z_marker);

        % Evenly spaced points in (0,1)
        u_x = linspace(0, 1, Nx + 2);  u_x = u_x(2:end-1);
        u_y = linspace(0, 1, Ny + 2);  u_y = u_y(2:end-1);
        u_z = linspace(0, 1, Nz + 2);  u_z = u_z(2:end-1);

        x_norm = norminv(u_x);
        y_norm = norminv(u_y);
        z_norm = norminv(u_z);

        % Rescale to volume bounds
        rescale = @(vals, lo, hi) lo + (vals - min(vals)) ./ ...
                                       max(eps, (max(vals) - min(vals))) * (hi - lo);

        x_conc = rescale(x_norm, volume(1,1), volume(1,2));
        y_conc = rescale(y_norm, volume(2,1), volume(2,2));
        z_conc = rescale(z_norm, volume(3,1), volume(3,2));

        [Xc, Yc, Zc] = meshgrid(x_conc, y_conc, z_conc);
        targetSpace  = [Xc(:), Yc(:), Zc(:)];

    otherwise
        error('generateTargetSpace:BadMode', ...
            'targetMode must be 1 (uniform) or 2 (normal). Got %g.', targetMode);
end
end
