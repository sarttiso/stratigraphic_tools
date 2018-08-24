% This function returns the strikes and dips corresponding to the normal 
% vectors describing planes of interest. Inverse of strdip2grad()
%
% IN:
% nx3 array, where n is the number of given strikes and dips. The
%   columns are the x,y,z components of the normal vector to the plane
%   described by the strike and dip.
%   x: increasing east
%   y: increasing north
%   z: increasing vertically
%
% OUT:
% strike: vector of strikes, measured with the right-hand rule
% dip: vectors of dips
%
% Adrian Tasistro-Hart, adrianraph-at-gmail.com, 24.08.2018

function [strike,dip] = grad2strdip(N)

x = N(:,1);
y = N(:,2);
z = N(:,3);

h = sqrt(x.^2+y.^2);

dip = 90 - atand(z./h);
strike = atand(x./y) - 90;

end