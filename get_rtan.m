%% Determine rpsi of the flux surface tangent to
% the beam trajectory.
%
% Input : jData -- JETPEAK data structure, obtained by
%                  running read_jData
%         clambda -- ratio of launching angle vs the one required
%                    to be tangent to the magnetic axis
%
% Output: rtan -- rpsi of the flux surface tangent to
%                 the beam trajectory
%
function rtan = get_rtan(jData, clambda)

lMag = get_lMag(jData);

Rtan = jData.Rmax * cos(lMag*clambda);

% Major radius of outboard midplane of each flux surface:
Rpsi = jData.Rmaj + jData.rpsi;

% Index of flux surface to which the beam is tangential:
iflxtan = find(abs(Rtan-Rpsi) == min(abs(Rtan-Rpsi)));

rtan = jData.rpsi(iflxtan);

% If tangential to the mangetic axis, then grids fail.
% Set manually:
if clambda == 1
    rtan = 0.0;
end

end
