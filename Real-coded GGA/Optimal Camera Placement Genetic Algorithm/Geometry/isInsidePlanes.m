function tf = isInsidePlanes(x, planes, tol, numVisible)
% Returns true if point x lies inside the intersection of all visible
% camera frustums. Each frustum is defined by a set of half-spaces
% (n_k . x >= d_k for inward-pointing normals, with a small tolerance tol
% which is expected to be a small negative number to admit boundary cases).
%
% The point is OUTSIDE camera i's frustum if ANY of its planes is violated
% (i.e. n_k . x - d_k < tol). The vectorised check below mirrors the
% intent of the commented-out scalar version: return false on the first
% violated plane, in any visible camera.
    for i= 1:numVisible
        nVecs = planes{i}(1:3,:); %normals stacked as columns
        dVals = planes{i}(4, :); %row of d offsets

        dotProd = nVecs' * x(:); %column: n_k . x for each plane k
        if any(dotProd - dVals' < tol) %any violated plane => x is outside this frustum
            tf = false;
            return;
        end
    end
    tf = true;
end












% 
% 
% Non-vectorised code
% 
%     for cc = 1:numel(planes) %loop through number of pyramids (cams) 
%         ps = planes{cc}; 
%         for k = 1:numel(ps) %loop through the plane surfaces 
%             if dot(ps(k).n, x) - ps(k).d < tol %if less than twiddle factor then not a vertex
%                 tf = false;
%                 return
%             end
%         end
%     end
%     tf = true;
% end