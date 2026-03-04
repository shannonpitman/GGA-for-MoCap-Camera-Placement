function targetSpace = generateTargetSpace(volume, targetMode, spacing)
% generates a Nx3 matrix of [x,y,z] target points inside the workspace
% volume

% Generate a uniform grid
x_marker = volume(1,1):spacing:volume(1,2);
y_marker = volume(2,1):spacing:volume(2,2);
z_marker = volume(3,1):spacing:volume(3,2);

switch targetMode
    case 1
        % Uniform Grid
        [X,Y,Z] = meshgrid(x_marker, y_marker, z_marker);
        targetSpace = [X(:),Y(:),Z(:)];
    case 2
        % Geneate normally-concentrated grid via inverse CDF warping
        Nx = numel(x_marker); % Number of points along each axis matches the uniform grid
        Ny = numel(y_marker);
        Nz = numel(z_marker);
    
        % Evenly spaced points in (0,1)
        u_x = linspace(0, 1, Nx + 2);  
        u_x = u_x(2:end-1);
        u_y = linspace(0, 1, Ny + 2);  
        u_y = u_y(2:end-1);
        u_z = linspace(0, 1, Nz + 2);  
        u_z = u_z(2:end-1);
    
        x_norm = norminv(u_x);  % inverse normal CDF 
        y_norm = norminv(u_y);
        z_norm = norminv(u_z);
    
        % Rescale to volume bounds: map [min(norm), max(norm)] -> [bound_lo, bound_hi]
        rescale = @(vals, lo, hi) lo + (vals - min(vals)) / (max(vals) - min(vals)) * (hi - lo);
    
        x_conc = rescale(x_norm, volume(1,1), volume(1,2));
        y_conc = rescale(y_norm, volume(2,1), volume(2,2));
        z_conc = rescale(z_norm, volume(3,1), volume(3,2));
    
        [Xc, Yc, Zc] = meshgrid(x_conc, y_conc, z_conc);
        targetSpace = [Xc(:), Yc(:), Zc(:)];
end