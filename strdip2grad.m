% This function returns the normal vector to a plane defined by the given
% strike and dip. Each component of the normal N corresponds to the x,y,z
% components of the gradient.
% x: increasing east
% y: increasing north
% z: increasing vertically

function N = strdip2grad(strike,dip)
    n = length(strike);
    strike = strike(:);
    dip = dip(:);

    d = [cosd(dip).*cosd(strike),...
         -cosd(dip).*sind(strike),...
         -sind(dip)];
    s = [sind(strike),...
         cosd(strike),...
         zeros(n,1)];
     
     N = cross(-s,d,2);
end