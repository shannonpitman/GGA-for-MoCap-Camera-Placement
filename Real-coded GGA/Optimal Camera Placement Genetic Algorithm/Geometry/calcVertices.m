function vertices = calcVertices(numVisible, visibleIdx, planes) 
% Vertices of the intersection of all visible camera frustums

    tol = -1e-9; % inside tolerance for boundary admission
    detThr = 1e-9; % "non-parallel" threshold for the 3-plane triple

    if numVisible < 2
        vertices = zeros(0, 3);
        return;
    end

    % Concatenate all visible cameras' planes into Nall (3xP), Dall (1xP)
    nPerCam = zeros(numVisible, 1);
    for i = 1:numVisible
        nPerCam(i) = size(planes{visibleIdx(i)}, 2);
    end
    P    = sum(nPerCam);
    Nall = zeros(3, P);
    Dall = zeros(1, P);
    col  = 0;
    for i = 1:numVisible
        c = nPerCam(i);
        Nall(:, col+1:col+c) = planes{visibleIdx(i)}(1:3, :);
        Dall(col+1:col+c) = planes{visibleIdx(i)}(4, :);
        col = col + c;
    end

    %Enumerate all triples of distinct planes. The triple INDEX set depends
    % only on P (number of planes = 4*numVisible), not on the geometry, and P
    % takes only a few values. Cache per P so nchoosek runs once per distinct
    % P per worker instead of on every point -- nchoosek/combs was the single
    % biggest cost in resUncertainty profiling (~160 s of self time).
    persistent tripleCache
    if isempty(tripleCache)
        tripleCache = {};
    end
    if numel(tripleCache) < P || isempty(tripleCache{P})
        tripleCache{P} = nchoosek(1:P, 3);
    end
    triples  = tripleCache{P};
    nTriples = size(triples, 1);

    % Three normals and offsets per triple
    aIdx = triples(:,1).'; 
    bIdx = triples(:,2).'; 
    cIdx = triples(:,3).';
    n_a = Nall(:, aIdx);  
    n_b = Nall(:, bIdx);  
    n_c = Nall(:, cIdx);
    d_a = Dall(aIdx);     
    d_b = Dall(bIdx);     
    d_c = Dall(cIdx);

    %Vectorised determinant filter
    nax = n_a(1,:); 
    nay = n_a(2,:); 
    naz = n_a(3,:);
    nbx = n_b(1,:); 
    nby = n_b(2,:); 
    nbz = n_b(3,:);
    ncx = n_c(1,:); 
    ncy = n_c(2,:); 
    ncz = n_c(3,:);

    bcyz = nby.*ncz - nbz.*ncy;
    bcxz = nbx.*ncz - nbz.*ncx;
    bcxy = nbx.*ncy - nby.*ncx;

    detA = nax.*bcyz - nay.*bcxz + naz.*bcxy;

    keep = abs(detA) > detThr;
    if ~any(keep)
        vertices = zeros(0, 3);
        return;
    end
    nax = nax(keep); 
    nay = nay(keep); 
    naz = naz(keep);
    nbx = nbx(keep); 
    nby = nby(keep); 
    nbz = nbz(keep);
    ncx = ncx(keep); 
    ncy = ncy(keep); 
    ncz = ncz(keep);
    d_a = d_a(keep); 
    d_b = d_b(keep); 
    d_c = d_c(keep);
    detA = detA(keep);
    bcyz = bcyz(keep); 
    bcxz = bcxz(keep); 
    bcxy = bcxy(keep);

    % Cramer's rule for x = A \ b across all kept triples
    detX = d_a.*bcyz -nay.*(d_b.*ncz - nbz.*d_c)+ naz.*(d_b.*ncy - nby.*d_c);

    detY = nax.*(d_b.*ncz - nbz.*d_c)- d_a.*bcxz+ naz.*(nbx.*d_c - d_b.*ncx);

    detZ = nax.*(nby.*d_c - d_b.*ncy)- nay.*(nbx.*d_c - d_b.*ncx)+ d_a.*bcxy;

    X = [detX; detY; detZ] ./ detA;

    % Iside-every-frustum test:(p, t) = n_p . x_t - d_p -> point t is outside if any p violates.
    violations = Nall.' * X - Dall.';
    inside     = ~any(violations < tol, 1);
    vertices = X(:, inside).';
end
