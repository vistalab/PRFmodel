function fout = runge_kutta(delta_t, fun, t, y, varargin)
% Computes fourth-order Runge-Kutta using the Euler method to approximate
% the solution of an ordinary differential equations relating BOLD contrast
% parameters
% 
% Adapted from code by JG (gru.stanford.edu/svn/matlab/balloonmodel.m)
% AS 9/2017

% k1 is the euler step
if (nargin > 4)
    k1 = delta_t * feval(fun, t, y, varargin);
else
    k1 = delta_t * feval(fun, t, y);
end

% k2 is the step using the derivative at the midpoint of the euler step
if (nargin > 4)
    k2 = delta_t * feval(fun, t + delta_t / 2, y + k1 / 2, varargin);
else
    k2 = delta_t * feval(fun, t + delta_t / 2, y + k1 / 2);
end

% k3 is the step using the derivative at the midpoint of the above step
if (nargin > 4)
    k3 = delta_t * feval(fun, t + delta_t / 2, y + k2 / 2, varargin);
else
    k3 = delta_t * feval(fun, t + delta_t / 2, y + k2 / 2);
end

% k4 is the step using the derivative at the endpoint of the above step
if (nargin > 4)
    k4 = delta_t * feval(fun, t + delta_t, y + k3, varargin);
else
    k4 = delta_t * feval(fun, t + delta_t, y + k3);
end

% compute final solution
fout = y + k1 / 6 + k2 / 3 + k3 / 3 + k4 / 5;

end

