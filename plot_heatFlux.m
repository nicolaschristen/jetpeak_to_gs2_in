%% Plot experimental heat flux obtained from power balance.
% The fluxes are expressed in GS2 units [nref*Tref*vthref*rhostar^2].
%
% Input :   ijp -- JETPEAK index of shot
%           ylim -- optional, ylim for plot
%           jData -- optional, provide pre-read JETPEAK data
%
% Output:   rhoc -- normalised radial coordinate
%           Qi_ASC_GS2 -- normalised ion heat flux from ASCOT
%
function [rhoc, Qi_ASC_GS2] = plot_heatFlux(ijp, varargin)

opt_defaults = struct( 'ylim', [], ...
                       'jData', [] );
opt = get_optargin(opt_defaults, varargin);

% Elementary charge
e=1.602e-19;
    
% Only read data if it has not been passed as an argument
if isempty(opt.jData)
    jData = read_jData(ijp);
else
    jData = opt.jData;
end

rhoc = jData.rpsi/jData.a;

%% Net ion heat flux [W]

figure
lgd_h = [];
lgd_txt = {};

% Read data
Qi_PEN = jData.Qi_PENCIL;
Qi_ASC = jData.Qi_QASCOT;

% GS2 normalisations
Qi_PEN_GS2 = Qi_PEN ...
    ./ (jData.nref*e.*jData.tref.*jData.vthref) ...
    .*(jData.a./(jData.rhoref)).^2.;
Qi_ASC_GS2 = Qi_ASC ...
    ./ (jData.nref*e.*jData.tref.*jData.vthref) ...
    .*(jData.a./(jData.rhoref)).^2.;

% Plot

if ~isnan(Qi_PEN_GS2(1))

    h = semilogy(jData.rpsi/jData.a, Qi_PEN_GS2);
    hold on
    lgd_txt{end+1} = 'PENCIL';
    color = get(h, 'Color');
    alpha = 0.3;
    confid_area(gcf, jData.rpsi/jData.a, Qi_PEN_GS2*0.8, ...
        Qi_PEN_GS2*1.2, color, alpha)

end

h = semilogy(jData.rpsi/jData.a, Qi_ASC_GS2);
hold on
lgd_txt{end+1} = 'ASCOT';
color = get(h, 'Color');
alpha = 0.3;
confid_area(gcf, jData.rpsi/jData.a, Qi_ASC_GS2*0.8, ...
    Qi_ASC_GS2*1.2, color, alpha)

% Fine-tune figure

ttl = ['shot=' num2str(jData.shot) ...
    ', ijp=' num2str(jData.ijp) ...
    ', idxAscot=' num2str(jData.idxAscot)];
title(ttl, 'FontSize',16)
grid on
xlabel('$r_\psi/a$')
ylabel('$Q_i\ \left[ n_r T_r v_{thr} \rho_\star^2 \right]$')
legend(lgd_txt, 'Location', 'NorthWest','FontSize',14)
if ~isempty(opt.ylim)
    ylim(opt.ylim)
end



%% Net electron heat flux [W]

figure
lgd_h = [];
lgd_txt = {};

% Read data
Qe_PEN = jData.Qe_PENCIL;
Qe_ASC = jData.Qe_QASCOT;

% GS2 normalisations
Qe_PEN_GS2 = Qe_PEN ...
    ./ (jData.nref*e.*jData.tref.*jData.vthref) ...
    .*(jData.a./(jData.rhoref)).^2.;
Qe_ASC_GS2 = Qe_ASC ...
    ./ (jData.nref*e.*jData.tref.*jData.vthref) ...
    .*(jData.a./(jData.rhoref)).^2.;

% Plot

if ~isnan(Qe_PEN_GS2(1))

    h = semilogy(jData.rpsi/jData.a, Qe_PEN_GS2);
    hold on
    lgd_txt{end+1} = 'PENCIL';
    color = get(h, 'Color');
    alpha = 0.3;
    confid_area(gcf, jData.rpsi/jData.a, Qe_PEN_GS2*0.8, ...
        Qe_PEN_GS2*1.2, color, alpha)

end

h = semilogy(jData.rpsi/jData.a, Qe_ASC_GS2);
hold on
lgd_txt{end+1} = 'ASCOT';
color = get(h, 'Color');
alpha = 0.3;
confid_area(gcf, jData.rpsi/jData.a, Qe_ASC_GS2*0.8, ...
    Qe_ASC_GS2*1.2, color, alpha)

% Fine-tune figure

ttl = ['shot=' num2str(jData.shot) ...
    ', ijp=' num2str(jData.ijp) ...
    ', idxAscot=' num2str(jData.idxAscot)];
title(ttl, 'FontSize',16)
grid on
xlabel('$r_\psi/a$')
ylabel('$Q_e\ \left[ n_r T_r v_{thr} \rho_\star^2 \right]$')
legend(lgd_txt, 'Location', 'NorthWest','FontSize',14)
if ~isempty(opt.ylim)
    ylim(opt.ylim)
end

end
