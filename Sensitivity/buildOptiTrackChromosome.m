function [chrom, calibID] = buildOptiTrackChromosome()
%BUILDOPTITRACKCHROMOSOME  Build the 1x42 chromosome for the OptiTrack rig.
%   Reproduces the lab measurements supplied for the 7-camera OptiTrack
%   ad-hoc rig and converts them into the GA's chromosome representation
%   (six [x y z rx ry rz] genes per camera: position + rotation vector).
%
%   The OptiTrack frame is Y-up and uses the OpenGL camera convention.
%   We transform to MATLAB's Z-up frame via T_transform, then post-multiply
%   each rotation by a 180-degree X rotation so the camera looks down +Z
%   in the convention used by setupCameras().
%
%   Source measurements: lab data exported from OptiTrack Motive,
%   reproduced verbatim from the user's earlier comparison script.

    T_transform = [1  0  0;
                   0  0 -1;
                   0  1  0];
    rotX180 = [1, 0, 0;
               0, cos(pi), -sin(pi);
               0, sin(pi),  cos(pi)];

    camPos = [
        -4.31837   3.16485  -1.40049;
        -4.30580   3.34875  -4.15416;
         4.41376   3.32004  -4.03407;
         4.38502   3.27992  -1.30325;
         4.36588   3.21471   1.41838;
         4.33365   3.33527   4.16674;
        -4.39346   3.30429   4.03519];

    R = cell(7,1);
    R{1} = [-0.0964006   0.553382  -0.82733;
            -0.0128461   0.830441   0.556959;
             0.99526     0.0643192 -0.0729462];
    R{2} = [-0.463454    0.348497  -0.814715;
             0.107535    0.934741   0.338667;
             0.879572    0.0693455 -0.470685];
    R{3} = [-0.383207   -0.381934   0.840999;
            -0.0928076   0.921818   0.376349;
            -0.918988    0.0661683 -0.388693];
    R{4} = [-0.0468672  -0.55851    0.828173;
            -0.0167187   0.829406   0.558395;
            -0.998761    0.0123245 -0.0482095];
    R{5} = [ 0.143307   -0.376961   0.915076;
             0.0465306   0.926163   0.374242;
            -0.988584   -0.0110526  0.150266];
    R{6} = [ 0.46651    -0.359219   0.808289;
             0.101597    0.929534   0.354465;
            -0.878662   -0.0832416  0.470132];
    R{7} = [ 0.770965    0.395039  -0.499556;
            -0.0314617   0.80705    0.589644;
             0.6361     -0.438878   0.634636];

    % Lens allocation. setupCameras() binds lens type to chromosome slot:
    % odd slots get the narrow (5.5 mm) lens, even slots the wide (3.5 mm).
    % The physical rig does NOT follow calibration order that way: the
    % four corner cameras (Motive IDs 2, 3, 6, 7) are narrow and the three
    % mid-wall cameras (IDs 1, 4, 5) are wide. So the calibration cameras
    % are re-ordered into slots such that every narrow camera lands in an
    % odd slot and every wide camera in an even slot. calibID(i) is the
    % Motive camera ID held in chromosome slot i (use it for labels).
    calibWide = [1 4 5];
    calibID   = [2 1 3 4 6 5 7];
    slotIsWide = mod(1:7, 2) == 0;
    assert(isequal(ismember(calibID, calibWide), slotIsWide), ...
        'buildOptiTrackChromosome:lensMismatch', ...
        'Slot order does not put the wide cameras (IDs 1, 4, 5) in even slots.');

    chrom = zeros(1, 6*7);
    for i = 1:7
        c = calibID(i);
        p_world = T_transform * camPos(c,:)';
        R_world = T_transform * R{c} * rotX180;
        chrom((i-1)*6 + (1:3)) = p_world';
        chrom((i-1)*6 + (4:6)) = rotmToGenes(R_world);   % rotation vector
    end
end
