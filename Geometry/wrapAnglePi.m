function y = wrapAnglePi(x)
%WRAPANGLEPI  Wrap angles (radians) into [-pi, pi).
%
%   y = wrapAnglePi(x) maps every element of x into the half-open interval
%   [-pi, pi). Used wherever a chromosome's Euler genes are edited so that
%   the result stays inside the search bounds used by runCameraOptimiser
%   (alpha, gamma in [-pi pi]).
%
%   Local reimplementation of wrapToPi so the project does not depend on
%   the Mapping / Robotics System toolboxes, which are not required
%   anywhere else in this codebase.

    y = mod(x + pi, 2*pi) - pi;
end
