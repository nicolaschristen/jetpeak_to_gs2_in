%% Script to reconstruct omega profiles for shots in the JETPEAK database.
%
%
% vvvvvvvvvvvvvvvvvvvvvvv %
% vvv USER PARAMETERS vvv %
% vvvvvvvvvvvvvvvvvvvvvvv %

% Paths to code
addpath('~/codes/matlab/tools/')
addpath('~/codes/jetpeak_v_gs2/')

% Shot index in JETPEAK
ijp = 950;

% Folder where GS2 fluxes and transport coeffs are saved
dataFolder = '~/data/gs2/prof_reconstruct/ijp_950/';
% Folder to save plots
plotFolder = [dataFolder 'userDepo_profiles/'];
% File containing GS2 fluxes
fluxFile = [dataFolder 'fluxes_gs2.csv'];
% File containing transport coeffs
tCoeffFile = [dataFolder 'transpCoeff.csv'];

% Use Pi/Q to reconstruct omega, instead of only Pi ?
use_Pi_over_Q = 0;
% flux ~ grad(psi)/|grad(Psi)| ? Else, flux ~ grad(x).
trinity_norm = 1;
% Aply GS2 normalisation to plotted quantites ?
nrm_gs2 = 0;
% Plot profiles of deposition and fluxes ?
plotverbose = 1;

% y-axis plotting limits
ylim_Q = [];
ylim_PI = [];
ylim_PI_over_Q = [];

% ^^^^^^^^^^^^^^^^^^^^^^^ %
% ^^^ USER PARAMETERS ^^^ %
% ^^^^^^^^^^^^^^^^^^^^^^^ %


format long

if nrm_gs2
    xlab = '$r_\psi/a$';
    ylab = '$\Omega_\phi$ [$v_{thr}$/$a$]';
else
    xlab = '$r_\psi$ [m]';
    ylab = '$\Omega_\phi$ [rad/s]';
end


%    ------------    %


%% First get the experimental profile from the JETPEAK database.


jData = read_jData(ijp, 'trinity_norm', trinity_norm);
r_expt = jData.rpsi;
om_expt = jData.omega;
if nrm_gs2
    r_expt = r_expt/jData.a;
    om_expt = om_expt*jData.a./jData.vthref;
end

%    ------------    %


%% Then compute the reconstructed profile.
% Process involves computing the transport coefficients
% and writing them to a file.


omegaFile = [dataFolder 'omegaReconstruct.csv'];
[r_sim, om_rcst, ~, ~] = omega_reconstruct( ijp, fluxFile, ...
    'fname_omega', omegaFile, ...
    'use_Pi_over_Q', use_Pi_over_Q, ...
    'nrm_gs2', nrm_gs2, ...
    'depoParams', [], ...
    'jData', jData, ...
    'trinity_norm', trinity_norm );
% Choose whether to normalise to GS2 units
if nrm_gs2
    r_sim = r_sim/jData.a;
    vthref_sim = interpol(jData.rpsi, jData.vthref, r_sim);
    om_rcst = om_rcst*jData.a./vthref_sim;
end


%    ------------    %


%% Parameters used to manually fit experimental deposition profiles
% with a Gaussian function.

depoParams.orig.SQi.width = 0.3;
depoParams.orig.SQi.nrm = 3.25e5;
depoParams.orig.SQi.mpos = 0.0;
depoParams.orig.SQi.skew = 0.0;
depoParams.orig.SQe.width = 0.75;
depoParams.orig.SQe.nrm = 3.75e5;
depoParams.orig.SQe.mpos = 0.0;
depoParams.orig.SQe.skew = 0.0;
depoParams.orig.SPi.width = 0.5;
depoParams.orig.SPi.nrm = 0.6;
depoParams.orig.SPi.mpos = 0.0;
depoParams.orig.SPi.skew = 0.0;



%    ------------    %


%% Compute omega with user-specified deposition profiles



% Cell-array with each element corresponding to a structure of user-specified
% parameters used to modify the deposition profiles.
usrParams = {};
% Cell-array of deposition profiles associated with usrParams
usrProfs = {};



% User-specified case #1:

% Sanity check. User-specified depositions are set to be equal to the fitted
% experimental values. If the GS2 simulations match the experiment well enough,
% this should give the same result as 'reconstructed' above.

% Same power as manual fit to experiment
depoParams.usr.cP = 1;
% Same launch angle as experiment
depoParams.usr.cl = 1;
% Same energy as experiment
depoParams.usr.cE = 0.5;

% Add these parameters to usrParams
usrParams{end+1} = depoParams.usr;

omegaFile = [dataFolder 'omega_doublePower.csv'];
[~, om_fittedDepo, srcQi, srcPI] = omega_reconstruct( ijp, tCoeffFile, ...
    'fname_omega', omegaFile, ...
    'use_Pi_over_Q', use_Pi_over_Q, ...
    'nrm_gs2', nrm_gs2, ...
    'depoParams', depoParams, ...
    'jData', jData, ...
    'trinity_norm', trinity_norm );

% Add corresponding deposition profiles to usrProfs
usrProfs{end+1}.srcQi = srcQi;
usrProfs{end}.srcPI = srcPI;

% Choose whether to normalise to GS2 units
if nrm_gs2
    om_fittedDepo = om_fittedDepo*jData.a./vthref_sim;
end

%    ------------    %


%% Plotting

% Plot omega profiles

figure
lgd_h = [plot(r_expt, abs(om_expt), 'LineWidth', 2)];
lgd_txt = {'Experiment'};
hold on
lgd_h(end+1) = plot(r_sim, abs(om_rcst), 'LineWidth', 2, 'Marker', 'o', ...
    'LineStyle', 'none', 'MarkerSize', 10);
lgd_txt{end+1} = 'GS2 reconstruction';
hold on
lgd_h(end+1) = plot(r_sim, abs(om_fittedDepo), 'LineWidth', 2, 'Marker', '+', ...
    'LineStyle', 'none', 'MarkerSize', 8);
lgd_txt{end+1} = 'GS2 + flux from experiment';
xlabel(xlab)
ylabel(ylab)
legend(lgd_h, lgd_txt, 'FontSize', 14, 'Location', 'SouthWest')
if use_Pi_over_Q
    ttl_omega = 'Method: $\chi_\varepsilon\Pi/\chi_LQ_i$';
else
    ttl_omega = 'Method: $\Pi/\chi_L$';
end
title(ttl_omega)
grid on

if plotverbose


    % Plot turbulent fluxes from the experiment and from GS2
    
    plot_fluxes( ijp, 'jData', jData, 'ylim_heat', ylim_Q, 'ylim_mom', ylim_PI, ...
                 'ylim_ratio', ylim_PI_over_Q, 'nrm_gs2', nrm_gs2, ...
                 'gs2_fluxFile', fluxFile, 'showTitle', 0, ...
                 'trinity_norm', trinity_norm );
            
    % Plot deposition profiles used

    plot_depo( ijp, 'jData', jData, 'origParams', depoParams.orig, ...
               'nrm_gs2', nrm_gs2 );

end
