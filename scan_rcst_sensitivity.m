%% Script to scan the sensitivity of the reconstruction process to
% the level of experimental fluxes given by ASCOT

ijp = 950;

dataFolder = '~/data/gs2/prof_reconstruct/ijp_950/scan_rcst_sensitivity/';
fluxFile = [dataFolder 'fluxes_gs2.csv'];
omegaFile = [dataFolder 'omegaReconstruct.csv'];
tCoeffFile = [dataFolder 'transpCoeff.csv'];

use_Pi_over_Q = 0;
nrm_gs2 = 0;

a.val = linspace(0.5,2,7);
a.idxmid = 3;
b.val = linspace(0.5,2,7);
b.idxmid = 3;
c.val = linspace(0.8,1.2,5);
c.idxmid = 3;

xlab = '$r_\psi$ [m]';
ylab = '$\Omega_\phi$ [s$^{-1}$]';

jData = read_jData(ijp);

[rpsi_rcst, om_rcst, ~, ~] = omega_reconstruct( ijp, fluxFile, ...
    'fname_omega', omegaFile, ...
    'use_Pi_over_Q', use_Pi_over_Q, ...
    'nrm_gs2', nrm_gs2, ...
    'depoParams', [], ...
    'jData', jData );
nr = numel(rpsi_rcst);

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

% Scan in alpha

a.N = numel(a.val);
a.om = zeros(nr,a.N);
a.cm = RdBu(a.N,'idxmid',a.idxmid);
for ia = 1:a.N
    [~, a.om(:,ia), ~, ~] = omega_reconstruct( ijp, tCoeffFile, ...
        'fname_omega', omegaFile, ...
        'use_Pi_over_Q', use_Pi_over_Q, ...
        'nrm_gs2', nrm_gs2, ...
        'depoParams', depoParams, ...
        'jData', jData, ...
        'fac_int', a.val(ia) );
end

figure
lgd_h = [plot(jData.rpsi, abs(jData.omega), 'k--')];
lgd_txt = {'Experiment'};
hold on
for ia = 1:a.N
    plot(rpsi_rcst, a.om(:,ia), 'Color', a.cm(ia,:))
    hold on
end
grid on
xlabel(xlab)
ylabel(ylab)
legend(lgd_h,lgd_txt)
colormap(a.cm);
cb = colorbar;
cb.Ticks = 0.5:0.5:2.0;
cb.FontSize = 14;
cb.TickLength = 0.0;
cb.Label.String = '$\alpha$';
cb.Label.FontSize = 24;
cb.Label.Interpreter = 'Latex';
caxis([0.5 2.0])

%% Scan in beta
%
%b.N = numel(b.val);
%b.om = zeros(nr,b.N);
%b.cm = RdBu(b.N,'idxmid',b.idxmid);
%for ib = 1:b.N
%    [~, b.om(:,ib), ~, ~] = omega_reconstruct( ijp, fluxFile, ...
%        'fname_omega', omegaFile, ...
%        'use_Pi_over_Q', use_Pi_over_Q, ...
%        'nrm_gs2', nrm_gs2, ...
%        'depoParams', [], ...
%        'jData', jData, ...
%        'fac_exp', b.val(ib) );
%end
%
%figure
%lgd_h = [plot(jData.rpsi, abs(jData.omega), 'k--')];
%lgd_txt = {'Experiment'};
%hold on
%for ib = 1:b.N
%    plot(rpsi_rcst, b.om(:,ib), 'Color', b.cm(ib,:))
%    hold on
%end
%grid on
%xlabel(xlab)
%ylabel(ylab)
%legend(lgd_h,lgd_txt)
%colormap(b.cm);
%cb = colorbar;
%cb.Ticks = 0.5:0.5:2.0;
%cb.FontSize = 14;
%cb.TickLength = 0.0;
%cb.Label.String = '$\beta$';
%cb.Label.FontSize = 24;
%cb.Label.Interpreter = 'Latex';
%caxis([0.5 2.0])
%
%% Scan in gamma
%
%c.N = numel(c.val);
%c.om = zeros(nr,c.N);
%c.cm = RdBu(c.N,'idxmid',c.idxmid);
%for ic = 1:c.N
%    [~, c.om(:,ic), ~, ~] = omega_reconstruct( ijp, fluxFile, ...
%        'fname_omega', omegaFile, ...
%        'use_Pi_over_Q', use_Pi_over_Q, ...
%        'nrm_gs2', nrm_gs2, ...
%        'depoParams', [], ...
%        'jData', jData, ...
%        'fac_BC', c.val(ic) );
%end
%
%figure
%lgd_h = [plot(jData.rpsi, abs(jData.omega), 'k--')];
%lgd_txt = {'Experiment'};
%hold on
%for ic = 1:c.N
%    plot(rpsi_rcst, c.om(:,ic), 'Color', c.cm(ic,:))
%    hold on
%end
%grid on
%xlabel(xlab)
%ylabel(ylab)
%legend(lgd_h,lgd_txt)
%colormap(c.cm);
%cb = colorbar;
%cb.Ticks = 0.8:0.2:1.2;
%cb.FontSize = 14;
%cb.TickLength = 0.0;
%cb.Label.String = '$\gamma$';
%cb.Label.FontSize = 24;
%cb.Label.Interpreter = 'Latex';
%caxis([0.8 1.2])
