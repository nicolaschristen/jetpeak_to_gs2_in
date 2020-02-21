%% Plot experimental heat and momentum fluxes, as well as their ratio
% All quantities have been normalised following GS2.
%
% Input :   ijp -- JETPEAK index of the shot
%           ylim_heat -- optional, plotting limits for heat flux
%           ylim_mom -- optional, plotting limits for mom flux
%           ylim_ratio -- optional, plotting limits for Pi/Q
%           jData -- optional, pre-read JETPEAK data structure
%
% Output:   -
%
function plot_fluxes(ijp, varargin)

opt_defaults = struct( 'ylim_heat', [], ...
                       'ylim_mom', [], ...
                       'ylim_ratio', [], ...
                       'jData', [] );
opt = get_optargin(opt_defaults, varargin);

if isempty(opt.jData)
    opt.jData = read_jData(ijp);
end

[rhoc, Qi_ASC_GS2] = plot_heatFlux(ijp, 'ylim', opt.ylim_heat, 'jData', opt.jData);
[~, PI_ASC_GS2] = plot_momFlux(ijp, 'ylim', opt.ylim_mom, 'jData', opt.jData);
    
% Plot

figure
lgd_h = [plot(rhoc, -1*PI_ASC_GS2./Qi_ASC_GS2)];
lgd_txt = {'ASCOT'};

% Fine-tune figure

ttl = ['shot=' num2str(opt.jData.shot) ...
    ', ijp=' num2str(opt.jData.ijp) ...
    ', idxAscot=' num2str(opt.jData.idxAscot)];
title(ttl, 'FontSize',16)
grid on
xlabel('$r_\psi/a$')
ylabel('$\sum_s -\Pi_s / Q_i$')
legend(lgd_h, lgd_txt, 'Location', 'NorthEast','FontSize',14)
if ~isempty(opt.ylim_ratio)
    ylim(opt.ylim_ratio)
end

end
