%% Return atan with output angle within [0:2pi]
%
% Input :   x,y -- such that tan(ang)=y/x
%
% Output:   ang -- atan(y/x), 0 <= ... <= 2pi
%
function ang=atan_2pi(x,y)
ang=0;
if x >= 0. && y >= 0.
    ang=atan(y/x);
elseif x < 0. && y >= 0.
    ang=pi+atan(y/x);
elseif x <= 0. && y < 0.
    ang=pi+atan(y/x);
elseif x > 0. && y < 0.
    ang=2.*pi+atan(y/x);
end
end