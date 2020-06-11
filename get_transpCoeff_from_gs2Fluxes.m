%% Compute tansport coefficients, given the turbulent fluxes at several
% raddi, and assuming that the desired transport coefficient are only
% weakly dependent on rotation, rotation shear, temperature and
% density. The transport coefficients are then written to a file.
%
% Input :   ijp --  shot index in JETPEAK DB
%           fname -- name of csv file containing computed fluxes,
%                    see example_files/fluxes.csv for format. The transport
%                    coefficients are saved in the same location as fname,
%                    in a file called 'transportCoeff.csv'
%           jData -- [optional] data structure obtained from JETPEAK, if
%                    read_jData has already been called.
%           trinity_norm -- [optional,0] if true, gs2 flux dotted with gradPsi
%                           else dotted with grad(x).
%
% Ouput:    flx -- table containing fluxes read from fname
%           tCoeff -- table containing the transport coefficients
%
function [flx, tCoeff] = get_transpCoeff_from_gs2Fluxes(ijp, fname, varargin)


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

%    ------------    %

% Elementary charge
cst.e = 1.602176634e-19;

%    ------------    %

% Read gs2 fluxes from file
flx = read_gs2Fluxes(ijp, fname, 'jData', jData, ...
                                 'trinity_norm', opt.trinity_norm );

%    ------------    %

% Compute transport coefficients

% First go back to standard signs for omega and its gradient,
% assuming omega does not change sign along the radial direction.
sgn_omega = sign(jData.omega(flx.ir_jData));

% Momentum pinch [m^2/s]
momPinch = -1*flx.PI_noGexb_gs2 .* flx.PINorm ./ ...
    (jData.mref.*jData.nref(flx.ir_jData).*jData.Rmaj(flx.ir_jData) ...
        .* abs(jData.omega(flx.ir_jData)))';
% Momentum diffusivity [m^2/s]
momDif = -1*(flx.PI_gs2-flx.PI_noGexb_gs2) .* flx.PINorm ./ ...
    (jData.mref.*jData.nref(flx.ir_jData).*jData.Rmaj(flx.ir_jData).^2 .* ...
        sgn_omega.*jData.domega_drpsi(flx.ir_jData))';
% Ion heat diffusivity [m^2/s]
% TODO: make this expression more precise ?
heatDif = -1*flx.Qi_gs2 .* flx.QNorm ./ ...
    ( cst.e*jData.nref(flx.ir_jData).*jData.dti_drpsi(flx.ir_jData) )';

% Write transport coefficients to file
% in same folder as file containing fluxes.
rpsi = flx.rpsi;
tCoeff = table(rpsi, momPinch, momDif, heatDif);
[fpath,~,~] = fileparts(fname);
writetable(tCoeff, [fpath '/transpCoeff.csv']);

end
