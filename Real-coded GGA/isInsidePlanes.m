function tf = isInsidePlanes(x, planes, tol, numVisible)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    for i= 1:numVisible
        nVecs = planes{i}(1:3,:); %normals stacked
        dVals = planes{i}(4, :); 

        dotProd = nVecs' *x(:);
        if all(dotProd-dVals' < tol) %Returns false if less  than twiddle factor
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