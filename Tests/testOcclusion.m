function tests = testOcclusion
% Unit tests for the dynamic-occlusion metric against Rahimian & Kearney
% (2017), Sec. 3.4 and Eq. 1. Run from the repo root:
%   matlab -batch "runtests('Tests')"
    tests = functiontests(localfunctions);
end

function setupOnce(~)
    addProjectPaths();
end

%% Eq. 1 penalties
function testNoCameras(tc)
    verifyEqual(tc, calculatePointOcclusion([], zeros(3, 0), 40, 140), 720);
end

function testOneCamera(tc)
    verifyEqual(tc, calculatePointOcclusion(1, [1; 0; 0], 40, 140), 360);
end

%% Paper Fig. 8: two triangulable cameras with horizontal separation beta
% are occluded over 180 + beta degrees.
function testPaperTwoCameraCase(tc)
    for beta = [40 60 90 120 140]
        V = [1, cosd(beta); 0, sind(beta); 0, 0];
        Q = calculatePointOcclusion(1:2, V, 40, 140);
        verifyEqual(tc, Q, 180 + beta, 'AbsTol', 1e-9, ...
            sprintf('beta = %d deg', beta));
    end
end

function testNonTriangulablePairIsFullyOccluded(tc)
    V = [1, cosd(20); 0, sind(20); 0, 0];   % 20 deg < 40 deg
    verifyEqual(tc, calculatePointOcclusion(1:2, V, 40, 140), 360, 'AbsTol', 1e-9);
end

%% A camera straight above the point sits on the occluder axis and is
% never occluded, so only the second camera's 180 deg of front-side counts.
function testCameraOnOccluderAxis(tc)
    V = [0, cosd(30); 0, 0; 1, sind(30)];   % 60 deg apart in 3D
    verifyEqual(tc, calculatePointOcclusion(1:2, V, 40, 140), 180, 'AbsTol', 1e-9);
end

%% Column order of the view vectors must not change Q (index-mismatch bug).
function testOrderInvariance(tc)
    rng(7);
    for trial = 1:50
        V = randomViewVectors(randi([2 7]));
        p = randperm(size(V, 2));
        verifyEqual(tc, calculatePointOcclusion(1:size(V,2), V(:, p), 40, 140), ...
            calculatePointOcclusion(1:size(V,2), V, 40, 140), 'AbsTol', 1e-9);
    end
end

%% Agreement with a dense brute-force sweep of occluder orientations.
function testMatchesBruteForceSweep(tc)
    rng(1);
    step = 0.01;
    for trial = 1:200
        n = randi([2 7]);
        V = randomViewVectors(n);
        Q = calculatePointOcclusion(1:n, V, 40, 140);
        verifyEqual(tc, Q, bruteForceQ(V, 40, 140, step), ...
            'AbsTol', 2*n*step + 1e-6, sprintf('trial %d', trial));
    end
end

%% Helpers
function V = randomViewVectors(n)
    V = randn(3, n);
    V(3, :) = 0.5*abs(V(3, :));
    V = V ./ vecnorm(V);
end

function Q = bruteForceQ(V, minAngle, maxAngle, step)
    va = mod(atan2d(V(2, :), V(1, :)), 360);
    phi = (step/2 : step : 360).';
    front = mod(va - phi, 360);
    front = front > 0 & front < 180;                    % orientations x cameras

    pairs = nchoosek(1:size(V, 2), 2);
    ang = acosd(max(-1, min(1, sum(V(:, pairs(:,1)) .* V(:, pairs(:,2)), 1))));
    good = pairs(ang >= minAngle & ang <= maxAngle, :);

    visible = false(size(phi));
    for k = 1:size(good, 1)
        visible = visible | (front(:, good(k,1)) & front(:, good(k,2)));
    end
    Q = step * sum(~visible);
end
