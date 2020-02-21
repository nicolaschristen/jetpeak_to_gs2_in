function [flx, tCoeff] = get_transpCoeff_from_gs2Fluxes(ijp, fname, varargin)


% Read optional input arguments
options_default = struct( 'jData', [] );
opt = get_optargin(options_default, varargin);

%    ------------    %

% Read data for this shot from JETPEAK
if isempty(opt.jData)
    jData = read_jData(ijp);
else
    jData = opt.jData;
end

% Read the fluxes from file
% flux > 0 -> radially outward
flx = readtable(fname);

% Determine indices in jData corresponding to rpsi in flx
ir_jData = zeros(numel(flx.rhoc),1);
for ir_flx = 1:numel(flx.rhoc)
    r = jData.a * flx.rhoc(ir_flx);
    ir_jData(ir_flx) = find(abs(jData.rpsi-r) == min(abs(jData.rpsi-r)));
end

%    ------------    %

% Elementary charge
e=1.602e-19;

% Distinguish normalised quantities from dimensionful ones
rpsi = flx.rhoc * jData.a;
PINorm = (jData.nref(ir_jData).*jData.mref.*jData.vthref(ir_jData).^2)' ...
    ./ (jData.a./(jData.rhoref(ir_jData)))'.^2.;
QNorm = (jData.nref(ir_jData)*e.*jData.tref(ir_jData) ...
        .*jData.vthref(ir_jData))' ...
    ./ (jData.a./(jData.rhoref(ir_jData)))'.^2.;

to_add = table(rpsi,PINorm,QNorm);
flx = [flx to_add];

%    ------------    %

% Compute transport coefficients

% First go back to standard signs for omega and its gradient,
% assuming omega does not change sign along the radial direction.
sgn_omega = sign(jData.omega(ir_jData));

% Momentum pinch [m^2/s]
momPinch = -1*flx.PI_noGexb_gs2 .* flx.PINorm ./ ...
    (jData.mref.*jData.nref(ir_jData).*jData.Rmaj(ir_jData).*abs(jData.omega(ir_jData)))';
% Momentum diffusivity [m^2/s]
momDif = -1*(flx.PI_gs2-flx.PI_noGexb_gs2) .* flx.PINorm ./ ...
    (jData.mref.*jData.nref(ir_jData).*jData.Rmaj(ir_jData).^2 .* ...
        sgn_omega.*jData.domega_drpsi(ir_jData))';

% Write transport coefficients to file
% in same folder as file containing fluxes.
tCoeff = table(rpsi,momPinch,momDif);
[fpath,~,~] = fileparts(fname);
writetable(tCoeff, [fpath '/transpCoeff.csv']);

end
