%% Plot experimental momentum flux obtained from ASCOT.
% The fluxes are expressed in GS2 units [nref*mref*vthref^2*Lref*rhostar^2].
%
% Input :   ijp -- vector of JETPEAK indices of shot
%           my_ylim -- optional, ylim for plot
%
% Output:   -
%
function plot_momFlux(ijp, my_ylim)

if nargin < 2
    my_ylim = [];
end

% Elementary charge
e=1.602e-19;
    
for ishot = 1:numel(ijp)

    jData = read_jData(ijp(ishot));
    
    %% Net ion heat flux [W]

    figure
    lgd_h = [];
    lgd_txt = {};

    % Read data
    PI_PEN = jData.PI_PENCIL;
    PI_ASC = jData.PI_ASCOT;

    % GS2 normalisations
    PI_PEN_GS2 = PI_PEN ...
        ./ (jData.nref.*jData.mref.*jData.vthref.^2) ...
        .*(jData.a./(jData.rhoref)).^2.;
    PI_ASC_GS2 = PI_ASC ...
        ./ (jData.nref.*jData.mref.*jData.vthref.^2) ...
        .*(jData.a./(jData.rhoref)).^2.;
    
    % Plot

    if ~isnan(PI_PEN_GS2(1))

        lgd_h(end+1) = semilogy(jData.rpsi/jData.a, PI_PEN_GS2);
        hold on
        lgd_txt{end+1} = 'PENCIL';
        color = get(lgd_h(end), 'Color');
        alpha = 0.3;
        confid_area(gcf, jData.rpsi/jData.a, PI_PEN_GS2*0.8, ...
            PI_PEN_GS2*1.2, color, alpha)

    end

    lgd_h(end+1) = semilogy(jData.rpsi/jData.a, PI_ASC_GS2);
    hold on
    lgd_txt{end+1} = 'ASCOT';
    color = get(lgd_h(end), 'Color');
    alpha = 0.3;
    confid_area(gcf, jData.rpsi/jData.a, PI_ASC_GS2*0.8, ...
        PI_ASC_GS2*1.2, color, alpha)

    % Fine-tune figure

    ttl = ['shot=' num2str(jData.shot) ...
        ', ijp=' num2str(jData.ijp) ...
        ', idxAscot=' num2str(jData.idxAscot)];
    title(ttl, 'FontSize',16)
    grid on
    xlabel('$r_\psi/a$')
    ylabel('$\sum_s\Pi_s \ \left[ n_r m_r v_{thr}^2 \rho_\star^2 \right]$')
    legend(lgd_h, lgd_txt, 'Location', 'NorthEast','FontSize',14)
    if ~isempty(my_ylim)
        ylim(my_ylim)
    end

end % loop over shots

end
