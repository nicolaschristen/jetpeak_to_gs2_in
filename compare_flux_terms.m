%% Plotting terms appearing in the transport equations
% to determine their relative sizes.


ijp = 950;

dataFolder = '~/data/gs2/prof_reconstruct/ijp_950/';
fluxFile = [dataFolder 'fluxes_gs2.csv'];
omegaFile = [dataFolder 'omegaReconstruct.csv'];
tCoeffFile = [dataFolder 'transpCoeff.csv'];

use_Pi_over_Q = 0;
nrm_gs2 = 0;

jData = read_jData(ijp);
flx = read_gs2Fluxes(ijp, fluxFile, 'jData', jData);

%% Compare terms in momentum equation

% Term from the turbulent momentum flux
PIterm = flx.PI;
% Term from the particle flux, due to change of frame,
% where <R^2 * Gamma>_psi is approximated by
% <Gamma>_psi * Rmag.
Gterm = flx.Gamma * jData.mref.*abs(jData.omega(flx.ir_jData))'*jData.Rmag^2;

figure
plot(flx.rpsi,PIterm)
hold on
plot(flx.rpsi,Gterm)
grid on
xlabel('$r_\psi$ [m]')
xlim([0 1])
legend({'$\Pi$','$m_i\Omega_\phi R_{mag}^2 \Gamma_i/$'},'Location','NorthWest')
title('Fluxes in momentum eq.')

%% Compare terms in heat equation

% External heating to ions, without taking into account
% the transfer wit electrons
S_ext = jData.srcE_i_QASCOT + jData.srcE_ie_QASCOT;
% Energy transfer term with electrons
S_transf = jData.srcE_ie_QASCOT;
% Turbulent heating, where <R^2 * Gamma>_psi is approximated by
% <Gamma>_psi * Rmag.
domega_dpsi = interpol(jData.psiflu,jData.omega,jData.psiflu,1);
domega_dpsi = domega_dpsi(flx.ir_jData)';
S_turb = abs(domega_dpsi) ./ jData.dx_dpsi(flx.ir_jData)' .* ...
    (flx.PI + jData.mref*jData.omega(flx.ir_jData)'*jData.Rmag^2.*flx.Gamma);

figure
plot(jData.rpsi,S_ext)
hold on
plot(jData.rpsi,S_transf)
hold on
plot(flx.rpsi,S_turb)
grid on
xlabel('$r_\psi$ [m]')
xlim([0 1])
legend('$S_p$','$S_{\Delta T}$','$H$')
title('Sources in heat eq.')

