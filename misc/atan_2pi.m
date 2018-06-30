%% Return atan with output theta in [0:2pi]
% Input :   x,y -- such that tan(theta)=y/x
function ang=atan_2pi(x,y)
ang=0;
if sign(x)==1 && sign(y)==1
    ang=atan(y/x);
elseif sign(x)==-1 && sign(y)==1
    ang=pi+atan(y/x);
elseif sign(x)==1 && sign(y)==-1
    ang=2.*pi+atan(y/x);
elseif sign(x)==-1 && sign(y)==-1
    ang=pi+atan(y/x);
end
end