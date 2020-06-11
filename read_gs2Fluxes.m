%% Read GS2 turbulent fluxes from a file
%
% Input :   ijp --  shot index in JETPEAK DB
%           fname -- name of csv file containing computed fluxes,
%                    see example_files/fluxes.csv for format. The transport
%                    coefficients are saved in the same location as fname,
%                    in a file called 'transportCoeff.csv'
%           jData -- [kw, []] data structure obtained from JETPEAK, if
%                    read_jData has already been called.
%           trinity_norm -- [kw, 0] if true, gs2 flux dotted with gradPsi
%                           else dotted with grad(x).
%
% Ouput:    flx -- table containing dimensionful fluxes read from fname
%
function flx = read_gs2Fluxes(ijp, fname, varargin)

% Read optional input arguments
options_default = struct( 'jData', [], ...
                          'trinity_norm', 0 );
opt = get_optargin(options_default, varargin);

%    ------------    %

% Read data for this shot from JETPEAK
if isempty(opt.jData)
    jData = read_jData(ijp, 'trinity_norm', opt.trinity_norm);
else
    jData = opt.jData;
end

% Read the fluxes from file
% flux > 0 -> radially outward
gs2dat = readtable(fname);

%    ------------    %

% Elementary charge
e=1.602e-19;

% Make quantities dimensionful
rpsi = gs2dat.rhoc * jData.a;

% Normalisation factors at simulated radii
PINorm = interpol(jData.rpsi, jData.PINorm, rpsi);
QNorm = interpol(jData.rpsi, jData.QNorm, rpsi);
GammaNorm = interpol(jData.rpsi, jData.GammaNorm, rpsi);

PI = gs2dat.PI_gs2.*PINorm;
PI_noGexb = gs2dat.PI_noGexb_gs2.*PINorm;
PI_noMach = gs2dat.PI_noMach_gs2.*PINorm;

Qi = gs2dat.Qi_gs2.*QNorm;
Qi_noGexb = gs2dat.Qi_noGexb_gs2.*QNorm;
Qi_noMach = gs2dat.Qi_noMach_gs2.*QNorm;

Qe = gs2dat.Qe_gs2.*QNorm;
Qe_noGexb = gs2dat.Qe_noGexb_gs2.*QNorm;
Qe_noMach = gs2dat.Qe_noMach_gs2.*QNorm;

Gamma = gs2dat.Gamma_gs2.*GammaNorm;
Gamma_noGexb = gs2dat.Gamma_noGexb_gs2.*GammaNorm;
Gamma_noMach = gs2dat.Gamma_noMach_gs2.*GammaNorm;

%    ------------    %

% Create a table with dimensionful quantities
flx = table( rpsi, PINorm, QNorm, GammaNorm, ...
                PI, PI_noGexb, PI_noMach, ...
                Qi, Qi_noGexb, Qi_noMach, ...
                Qe, Qe_noGexb, Qe_noMach, ...
                Gamma, Gamma_noGexb, Gamma_noMach );

end
