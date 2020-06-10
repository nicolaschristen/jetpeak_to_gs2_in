%% Script to reconstruct omega profiles for shots in the JETPEAK database.
%
%
% vvvvvvvvvvvvvvvvvvvvvvv %
% vvv USER PARAMETERS vvv %
% vvvvvvvvvvvvvvvvvvvvvvv %

addpath('~/codes/jetpeak_v_gs2/')
format long

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
use_Pi_over_Q = 1;
% Pi ~ grad(x) ? Else, Pi ~ grad(psi)/|grad(Psi)|.
trinity_norm = 1;
% Aply GS2 normalisation to final quantites ?
nrm_gs2 = 0;
% Plot profiles of deposition and fluxes ?
plotverbose = 1;

% y-axis plotting limits
ylim_Q = [];
ylim_PI = [];
ylim_PI_over_Q = [];

if nrm_gs2
    xlab = '$r_\psi/a$';
    ylab = '$\Omega_\phi$ [$v_{thr}$/$a$]';
else
    xlab = '$r_\psi$ [m]';
    ylab = '$\Omega_\phi$ [rad/s]';
end

% ^^^^^^^^^^^^^^^^^^^^^^^ %
% ^^^ USER PARAMETERS ^^^ %
% ^^^^^^^^^^^^^^^^^^^^^^^ %



%    ------------    %


%% First get the experimental profile from the JETPEAK database.


jData = read_jData(ijp);
r_expt = jData.rpsi;
om_expt = jData.omega;
if nrm_gs2
    r_expt = r_expt/jData.a;
    om_expt = om_expt*jData.a./jData.vthref;
end

%    ------------    %


%% Then compute the reconstructed profile.
% At the same time, compute the transport coefficients
% and write them to a file.


omegaFile = [dataFolder 'omegaReconstruct.csv'];
[r, om_rcst, ~, ~] = omega_reconstruct( ijp, fluxFile, ...
    'fname_omega', omegaFile, ...
    'use_Pi_over_Q', use_Pi_over_Q, ...
    'nrm_gs2', nrm_gs2, ...
    'depoParams', [], ...
    'jData', jData, ...
    'trinity_norm', trinity_norm );
% Choose whether to normalise to GS2 units
if nrm_gs2
    r = r/jData.a;
    for ir = 1:numel(r)
        rtmp = r(ir);
        ir_jetpeak2rcst = find(abs(jData.rpsi-rtmp) == min(abs(jData.rpsi-rtmp)));
        om_rcst(ir) = om_rcst(ir)*jData.a/jData.vthref(ir_jetpeak2rcst);
    end
end

%    ------------    %


%% Compute omega with user-specified deposition profiles


% Cell-array containing structures with each of the user-set deposition
% parameters
usrParams = {};
% Cell-array with associated deposition profiles
usrProfs = {};

% Specify width and norm of manual fit for the
% experimental deposition profiles

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

% Sanity check: specify the deposition profile such that they fit the
% experimental profiles. Should give same result as 'reconstructed' above.

% Same power as manual fit to experiment
depoParams.usr.cP = 1;
% Same launch angle as experiment
depoParams.usr.cl = 1;
% Same energy as experiment
depoParams.usr.cE = 0.5;

% Add these parameters to the usrParams
usrParams{end+1} = depoParams.usr;

omegaFile = [dataFolder 'omega_doublePower.csv'];
[~, om_fittedDepo, srcQi, srcPI] = omega_reconstruct( ijp, tCoeffFile, ...
    'fname_omega', omegaFile, ...
    'use_Pi_over_Q', use_Pi_over_Q, ...
    'nrm_gs2', nrm_gs2, ...
    'depoParams', depoParams, ...
    'jData', jData, ...
    'trinity_norm', trinity_norm );

% Add user-specified deposition profiles to the usrProfs
usrProfs{end+1}.srcQi = srcQi;
usrProfs{end}.srcPI = srcPI;

% Choose whether to normalise to GS2 units
if nrm_gs2
    for ir = 1:numel(r)
        rtmp = r(ir);
        ir_jetpeak2fit = find(abs(jData.rpsi-rtmp) == min(abs(jData.rpsi-rtmp)));
        om_fittedDepo(ir) = om_fittedDepo(ir)*jData.a/jData.vthref(ir_jetpeak2fit);
    end
end

%    ------------    %


%% Plotting

% Plot omega profiles

figure
lgd_h = [plot(r_expt, abs(om_expt), 'LineWidth', 2)];
lgd_txt = {'Experiment'};
hold on
lgd_h(end+1) = plot(r, abs(om_rcst), 'LineWidth', 2, 'Marker', 'o', ...
    'LineStyle', 'none', 'MarkerSize', 10);
lgd_txt{end+1} = 'GS2 reconstruction';
hold on
lgd_h(end+1) = plot(r, abs(om_fittedDepo), 'LineWidth', 2, 'Marker', '+', ...
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
