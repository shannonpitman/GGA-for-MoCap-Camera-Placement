function tf = isInsidePlanes(x, planes, tol, numVisible)
% True if x lies inside the intersection of all visible camera frustums.
    if nargin < 4 || isempty(numVisible)
        numVisible = numel(planes);
    end

    % Normalise input to 3 x T column-points.
    if isvector(x)
        X = x(:);
    else
        X = x;
        if size(X, 1) ~= 3 && size(X, 2) == 3
            X = X.';
        end
    end

    % Build the concatenated plane matrix once.
    nPerCam = zeros(numVisible, 1);
    for i = 1:numVisible
        nPerCam(i) = size(planes{i}, 2);
    end
    P = sum(nPerCam);
    if P == 0
        tf = true(1, size(X, 2));
        if size(X, 2) == 1, tf = tf(1); end
        return;
    end

    Nall = zeros(3, P);
    Dall = zeros(1, P);
    col  = 0;
    for i = 1:numVisible
        c = nPerCam(i);
        Nall(:, col+1:col+c) = planes{i}(1:3, :);
        Dall(col+1:col+c) = planes{i}(4, :);
        col = col + c;
    end

    % violations(p, t) = n_p . x_t - d_p ; point t is outside if any p violates.
    violations = Nall.' * X - Dall.';  % P x T
    tf = ~any(violations < tol, 1); % 1 x T

    if size(X, 2) == 1
        tf = tf(1);
    end
end
