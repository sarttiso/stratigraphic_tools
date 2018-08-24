% This function performs a least squares fit of a plane to a set of points.
%
% IN:
% x: x-coordinates of points
% y: y-coordinates of points
% z: z-coordinates of points
%
% OUT:
% N: 1x3 vector with components of vector normal to the fitted plane
% phi: sum of squared misfits between z coordinates of points and plane

function [N,phi] = fitplane(x,y,z)

% demean
x = x-mean(x);
y = y-mean(y);
z = z-mean(z);

% our model is as follows:
% equation for plane: ax+by+cz+d = 0
% we assume we have noisy measurements of z which would otherwise lie on
% the plane, so we rearrange the plane equation to solve for z':
% z' = -(d+ax+by)/c
% and we rename a/c = alpha, b/c = beta, d/c = gamma
% then we set up the least squares system:
% want to minimize sum (z'-z)^2 over alpha, beta, and gamma
% which means differentiating against those parameters. when you work it
% out, the system becomes:
% sum_i alpha*xi^2 + beta*yi*xi + gamma*xi = -sum_i xi*zi
% sum_i alpha*xi*yi + beta*yi^2 + gamma*yi = -sum_i yi*zi
% sum_i alpha*xi + beta*yi + gamma = -sum_i zi
% which readily decomposes into the linear system below.
A = [sum(x.^2) sum(x.*y) sum(x);
     sum(x.*y) sum(y.^2) sum(y);
     sum(x)    sum(y)    1];    % model A, grad(sum of squared error)
d = -[sum(x.*z);
     sum(y.*z);
     sum(z)];

% solve linear system
m = inv(A'*A)*A'*d;
alpha = m(1);
beta = m(2);
gamma = m(3);

% then, we see that the normal vector to the plane is just [a,b,c] which is
% proportional to [alpha,beta,1]:
N = [alpha,beta,1];

% misfit is easily computed:
phi = sum( (z + gamma + alpha*x + beta*y).^2);

end