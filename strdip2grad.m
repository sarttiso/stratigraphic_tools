% This function returns the normal vector to a plane defined by the given
% strike and dip. Each component of the normal N corresponds to the x,y,z
% components of the gradient.
%
% IN:
% strike: vector of strikes, measured with the right-hand rule, positive
%   clockwise from north
% dip: vectors of dips, must have same length as strike
%
% OUT:
% N: nx3 array, where n is the number of given strikes and dips. The
%   columns are the x,y,z components of the normal vector to the plane
%   described by the strike and dip.
%   x: increasing east
%   y: increasing north
%   z: increasing vertically
%
% Adrian Tasistro-Hart, adrianraph-at-gmail.com, 24.08.2018

function N = strdip2grad(strike,dip)

n = length(strike);
% make columns
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