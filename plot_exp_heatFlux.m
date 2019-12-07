%% Plot experimental heat flux obtained from power balance.
% The fluxes are expressed in GS units [nref*Tref*vthref*rhostar^2].
%
% Input :   ijp -- vector of JETPEAK indices of shot
%           my_ylim -- optional, ylim for plot
%
% Output:   -
%
function plot_exp_heatFlux(ijp, my_ylim)

if nargin < 2
    my_ylim = [];
end

load('~/codes/jetpeak_and_gs2/databases/TRANSP_2017_3.mat','TRANSP')
load('~/codes/jetpeak_and_gs2/databases/JETPEAK_2019_10.mat','Q')

figure
lgd_h = [];
lgd_txt = {};

for ishot = 1:numel(ijp)

    jData = read_jData(ijp(ishot));
    
    % Elementary charge
    e=1.602e-19;
    
    % Flux surface areas from TRANSP [m^{-2}]
    A_psi = 1e-4*TRANSP.T.SURF(TRANSP.JPI==ijp(ishot),:);
    
    % Read heatfluxes
    % Net ion heat flux from PENCIL & PION [W]
    Qi_pencil = Q.NBPICH.Qi(ijp(ishot),:);
    Qi_pencil = interpol(jData.sqrt_psin_chain2,Qi_pencil,jData.sqrt_psin_TRANSP);
    
    % GS2 normalisations
    Qi_pencil_GS2 = Qi_pencil ...
        ./ (jData.nref*e.*jData.tref.*jData.vthref) ...
        .*(jData.a./(jData.rhoref)).^2./A_psi;
    
    h = semilogy(jData.rpsi/jData.a, Qi_pencil_GS2);
    lgd_h(end+1) = h;
    lgd_txt{end+1} = ['Shot $\#' num2str(jData.shot) '$'];
    hold on
    color = get(h, 'Color');
    alpha = 0.3;
    confid_area(gcf, jData.rpsi/jData.a, Qi_pencil_GS2*0.8, ...
        Qi_pencil_GS2*1.2, color, alpha)

end % loop over shots

grid on
xlabel('$r_\psi/a$')
ylabel('$Q_{i}\ \left[ n_r T_r v_{thr} \rho_\star^2 \right]$')
legend(lgd_h, lgd_txt, 'Location', 'NorthWest','FontSize',14)
if ~isempty(my_ylim)
    ylim(my_ylim)
end

end
