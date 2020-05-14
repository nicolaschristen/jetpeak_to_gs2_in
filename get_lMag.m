%% Determine the launchaing angle for the beam to be
% tangent to the magnetic axis.
%
% Input : jData -- JETPEAK data structure, obtained by
%                  running read_jData
%
% Output: lMag -- launching angle for the beam to be
%                 tangent to the magnetic axis
%
function lMag = get_lMag(jData)

lMag = acos(jData.Rmag/jData.Rmax);

end
