function E = optimiseEllipsoid(axes, transformed_V, axes0, weight) 
    a = axes(1);
    b = axes(2);
    c = axes(3);

    xV = transformed_V(1,:);
    yV = transformed_V(2,:);
    zV = transformed_V(3,:);

    Tau = (xV.^2)/(a^2) + (yV.^2)/(b^2) + (zV.^2)/(c^2);
    ray = 1 ./ sqrt(Tau); %zE/zV
    %compute intersection with ellipsoid with ray piercing
    xE = ray .* xV; 
    yE = ray .* yV;
    zE = ray .* zV;
    
    dist_squared = (xV- xE).^2 +(yV-yE).^2 +(zV-zE).^2; %if vertex is already inside (then distance ==0)
    E1 = sum(dist_squared); % expand axes to cover vertices 
    E2 = sum((abs(axes0)-axes).^2); %don't inflate one axis unnecessarily 
    E = E1 + weight*E2;
end
