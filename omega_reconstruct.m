%% Reconstruct a rotation profile, given the turbulent fluxes at several
% raddi, and assuming that the corresponding transport coefficient are
% only weakly dependent on rotation, rotation shear, temperature and
% density. The function can be used in two ways:
%
%     1. Turbulent fluxes have been computed in simulations, and the
%        transport coefficients are calculated based on those fluxes.
%
%     2. Turbulent fluxes are computed based on user-specified deposition
%        profiles, and the transport coeffients are from a previous but
%        similar case where method 1. was used.
%
% Input :   ijp --  shot index in JETPEAK DB
%           fname -- [usage 1] name of csv file containing computed fluxes,
%                              see example_files/fluxes.csv for format.
%                    [usage 2] name of csv file containing computed
%                              transport coefficients, see
%                              example_files/transpCoeff.csv for format.
%           use_Pi_over_Q -- [kw, 0] if true, Pi/Q is used in reconstruction.
%                            Else, Pi is used.
%           depoParams -- [kw, []] if empty, sources from the experiment are used.
%                         To specify depositions with a skewed gaussian:
%                             * See set_useDepo.m for more details.
%                             * store width, nrm, mpos & skew of the fit to
%                               experimental profiles as fields of depoParams.orig
%                             * store cP, cl & cE as fields of depoParams.usr
%           fname_omega -- [kw, 'omegaReconstruct.csv'] file name to save
%                          rotation profiles.
%           nrm_gs2 -- [kw, 0] normalise plotted quantities to GS2 units
%           jData -- [kw, []] structure containing pre-read JETPEAK data
%           plotverbose -- [kw, 0] plot fluxes and deposition profiles
%           fac_int -- [kw, 1] factor scaling Pi or Pi/Q (for debugging)
%           trinity_norm -- [optional,0] if true, gs2 flux dotted with gradPsi
%                           else dotted with grad(x).
%
% Ouput:    -
%
function [rpsi, omega, srcQi, srcPI] = omega_reconstruct(ijp, fname, varargin)


% Read optional input arguments
options_default = struct( 'use_Pi_over_Q', 0, ...
                          'depoParams', [], ...
                          'fname_omega', 'omegaReconstruct.csv', ...
                          'nrm_gs2', 0, ...
                          'jData', [], ...
                          'plotverbose', 0, ...
                          'fac_int', 1, ...
                          'trinity_norm', 0 );
opt = get_optargin(options_default, varargin);

% Flag to tell if the deposition profile is from the experiment, or
% if it is user-specified
if isempty(opt.depoParams)
    use_origDepo = 1;
else
    use_origDepo = 0;
end

% Elementary charge
cst.e = 1.602176634e-19;

%    ------------    %

% Read data for this shot from JETPEAK
if isempty(opt.jData)
    jData = read_jData(ijp);
else
    jData = opt.jData;
end

% Either read the fluxes from simulations,
% use them to compute the transport coefficients,
% and save those coefficients to a file.
if use_origDepo
    [flx, tCoeff] = get_transpCoeff_from_gs2Fluxes(ijp, fname, 'jData', jData, ...
                                                   'trinity_norm', opt.trinity_norm );
% Or compute the fluxes based on user-specified profiles,
% and read pre-computed transport coefficients from a file.
else
    % Get the user-specified deposition profiles
    [usrDepo, origDepo] = set_userDepo( ijp, opt.depoParams.orig, opt.depoParams.usr, ...
                                        'jData', jData, 'plotverbose', opt.plotverbose, ...
                                        'nrm_gs2', opt.nrm_gs2 );
    % Compute the assiocated fluxes
    flx.rpsi = usrDepo.rpsi;
    flx.PI = flux_from_source(jData.psiflu, jData.dV, jData.dx_dpsi, usrDepo.srcPi);
    flx.Qi = flux_from_source(jData.psiflu, jData.dV, jData.dx_dpsi, usrDepo.srcQi);

    % Read transport coefficients from file
    tCoeff = readtable(fname);
end

%    ------------    %

% Define an rpsi interval on which to reconstruct omega
rpsiIn = min(tCoeff.rpsi); % innermost point
rpsiOut = max(tCoeff.rpsi); % outermost point
% Index in jData corresponding to rpsiOut
iout_jData = find(abs(rpsiOut-jData.rpsi) == min(abs(rpsiOut-jData.rpsi)));
nr = 1000;
recon.rpsi = linspace(rpsiIn,rpsiOut,nr);

%    ------------    %

% For quantities determined from simulations, use linear interpolation.
recon.momPinch = interp1(tCoeff.rpsi,tCoeff.momPinch,recon.rpsi);
recon.momDif = interp1(tCoeff.rpsi,tCoeff.momDif,recon.rpsi);
recon.heatDif = interp1(tCoeff.rpsi,tCoeff.heatDif,recon.rpsi);
% For user-set or experimental quantities, use splines.
recon.nref = interpol(jData.rpsi,jData.nref,recon.rpsi);
recon.Rmaj = interpol(jData.rpsi,jData.Rmaj,recon.rpsi);
recon.dti_drpsi = interpol(jData.rpsi,jData.dti_drpsi,recon.rpsi);
% If fluxes have been computed from simulations, then interpolate linearly
if use_origDepo
    recon.PI = interp1(flx.rpsi,flx.PI_gs2.*flx.PINorm,recon.rpsi);
    recon.Qi = interp1(flx.rpsi,flx.Qi_gs2.*flx.QNorm,recon.rpsi);
% If fluxes are computed from user-specified deposition profiles,
% then use splines.
else
    recon.PI = interpol(flx.rpsi,flx.PI,recon.rpsi);
    recon.Qi = interpol(flx.rpsi,flx.Qi,recon.rpsi);
end

%    ------------    %

% Integrate inwards using trapezoidal rule

dr = recon.rpsi(2) - recon.rpsi(1);

% First compute exponent
g = zeros(1,nr);
if ~opt.use_Pi_over_Q
    intgrd = recon.momPinch./(recon.momDif.*recon.Rmaj);
else
    intgrd = recon.momPinch./(recon.momDif.*recon.Rmaj);
end
for ir = nr-1:-1:1
    g(ir) = g(ir+1) + (intgrd(ir)+intgrd(ir+1)) * dr/2.0;
end

% Then compute main integral
f = zeros(1,nr);
if ~opt.use_Pi_over_Q
    intgrd = recon.PI./(recon.momDif.*recon.nref.*jData.mref.*recon.Rmaj.^2) .* exp(-1*g);
else
    intgrd = cst.e*recon.dti_drpsi./(jData.mref.*recon.Rmaj.^2) ...
                 .* recon.heatDif./recon.momDif .* recon.PI./recon.Qi .* exp(-1*g);
end
for ir = nr-1:-1:1
    f(ir) = f(ir+1) + (intgrd(ir)+intgrd(ir+1)) * dr/2.0;
end
% Apply optional scaling factor to integral (default=1)
f = opt.fac_int * f;

%    ------------    %

% Get and save expression for omega(rpsi)

% Edge BC for the angular frequency
omegaOut = abs(jData.omega(iout_jData));

% Taking care of sign changes when using only Pi,
% or Pi/Qi:
if ~opt.use_Pi_over_Q
    s = 1;
else
    s = -1;
end

% Compute reconstructed omega
recon.omega = (omegaOut*ones(1,nr) + s*f) .* exp(g);

% Get nearest neighbour on recon grid to the rpsi grid from the simulations
nr_tCoeff = numel(tCoeff.rpsi);
omega = zeros(nr_tCoeff,1);
for ir = 1:nr_tCoeff
    r = tCoeff.rpsi(ir);
    ir_recon2tCoeff = find(abs(recon.rpsi-r) == min(abs(recon.rpsi-r)));
    omega(ir) = recon.omega(ir_recon2tCoeff);
end

% Write reconstructed rotation profile to file
% in same folder as the file specified by fname.
rpsi = tCoeff.rpsi;
omRecon = table(rpsi,omega);
writetable(omRecon, opt.fname_omega);

% Give source terms as output
if ~use_origDepo
    srcQi = usrDepo.srcQi;
    srcPI = usrDepo.srcPi;
else
    srcQi = [];
    srcPI = [];
end

%    ------------    %

% Plotting

if opt.plotverbose

    figure
    h = [plot(omRecon.rpsi,omRecon.omega,'.-','MarkerSize',20,'LineWidth',1)];
    lgd = {'reconstructed'};
    if isempty(opt.depoParams)
        hold on
        h(end+1) = plot(jData.rpsi,abs(jData.omega));
        lgd{end+1} = 'experiment';
    end
    legend(h,lgd)
    xlabel('$r_\psi / a$')
    ylabel('$\Omega_\phi$ [rad/s]')
    grid on

end

end
