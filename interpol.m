%% Interpolate and take derivative (at specified order) of function f
%
% Input :   x -- points at which function f is known
%           fval -- values of f at x
%           xout (= x) -- points at which the new function is evaluated
%           order (= 0) -- order of the derivative that should be applied to f
%
% Output:   fun_der -- array containing order-th derivative of f,
%                      evaluated at x_out
%
function fun_der = interpol(x, fval, x_out, order)
if nargin < 3
    x_out=x;
end
if nargin < 4
    order=0;
end
fun_pp=spline(x,fval); % get piecewise polynomial expression
fun_der_pp=fnder(fun_pp,order); % take the derivative of expression
fun_der=ppval(fun_der_pp,x_out); % evaluate derivative
end