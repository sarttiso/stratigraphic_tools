function [strike, dip, m, phi] = strikedip(x,y,z)

%% Least Squares Strikes and Dips from Bed Traces 
% Model elevations as resulting from f(x,y) -> (x,y,f(x,y)) so that the
% plane we're looking for is z = Ax + By + C. We want to minimize misfit
% then which is misfit = sum((Ax + By + C - z)^2), so we set the gradient
% equal to zero. Makes problem linear with d = A*m, m = [A B C], A is
% complicated, and d = [sum(xi*zi) sum(yi*zi) sum(zi)]
%%
x = x-mean(x);
y = y-mean(y);
z = z-mean(z);

A = [sum(x.^2) sum(x.*y) sum(x);
     sum(x.*y) sum(y.^2) sum(y);
     sum(x)    sum(y)    1];    % model A, grad(sum of squared error)
d = [sum(x.*z);
     sum(y.*z);
     sum(z)];
 
m = inv(A'*A)*A'*d;
A = m(1);
B = m(2);
C = m(3);

n = [A B -1];
n = n/norm(n);    % unit normal vector

% dip is [dz/dx=A dz/dy=B t] and orthogonal to unit normal so dot(d,n) = 0
% which means [A B -1]*[A B t] = 0 = A^2 + B^2 -t so t = A^2+B^2
d = [A B A^2+B^2];
d = d/norm(d);
dh = [A B 0];  % horizontal projection of dip
dip = rad2deg(atan2(norm(cross(d,dh)),dot(d,dh)));

% strike is the cross product between normal and dip
s = cross(n,d);
strike = cart2compass(s(1),s(2));

phi = sum((z-(A*x+B*y+C)).^2);

end