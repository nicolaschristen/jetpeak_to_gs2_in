%% Plot experimental heat and momentum fluxes, as well as their ratio
% All quantities have been normalised following GS2.
%
% Input :   ijp -- JETPEAK index of the shot
%           ylim_heat -- [kw, []] plotting limits for heat flux
%           ylim_mom -- [kw, []] plotting limits for mom flux
%           ylim_ratio -- [kw, []] plotting limits for Pi/Q
%           jData -- [kw, []] pre-read JETPEAK data structure
%           gs2_fluxFile -- [kw, ''] file containing GS2 fluxes
%                           to be plotted
%           nrm_gs2 -- [kw, 0] normalise plotted quantities to GS2 units
%           showTitle -- [kw, 1] add title to plots
%           showAllCodes -- [kw, 0] plot other deposition codes than
%                           ASCOT, eg PENCIL
%           trinity_norm -- [kw, 0] if true, gs2 flux dotted with gradPsi
%                           else dotted with grad(x).
%
% Output:   -
%
function plot_fluxes(ijp, varargin)

opt_defaults = struct( 'ylim_heat', [], ...
                       'ylim_mom', [], ...
                       'ylim_ratio', [], ...
                       'jData', [], ...
                       'gs2_fluxFile', '', ...
                       'nrm_gs2', 0 , ...
                       'showTitle',1, ...
                       'showAllCodes', 0, ...
                       'trinity_norm', 0 );
opt = get_optargin(opt_defaults, varargin);

if isempty(opt.jData)
    opt.jData = read_jData(ijp);
end

% Plot Q

[~, ~, ~] = plot_heatFlux(ijp, 'ylim', opt.ylim_heat, 'jData', opt.jData, ...
                                   'gs2_fluxFile', opt.gs2_fluxFile, ...
                                   'nrm_gs2', opt.nrm_gs2, ...
                                   'showTitle', opt.showTitle, ...
                                   'showAllCodes', opt.showAllCodes, ...
                                   'trinity_norm', opt.trinity_norm );

% Plot Pi

[~, ~] = plot_momFlux(ijp, 'ylim', opt.ylim_mom, 'jData', opt.jData, ...
                                   'gs2_fluxFile', opt.gs2_fluxFile, ...
                                   'nrm_gs2', opt.nrm_gs2, ...
                                   'showTitle', opt.showTitle, ...
                                   'showAllCodes', opt.showAllCodes, ...
                                   'trinity_norm', opt.trinity_norm );
    
% Read GS2 fluxes form file

if ~isempty(opt.gs2_fluxFile)
    flx = read_gs2Fluxes( ijp, opt.gs2_fluxFile, 'jData', opt.jData, ...
                          'trinity_norm', opt.trinity_norm );
end

% Plot Pi/Q

figure

if opt.nrm_gs2
    xvar = opt.jData.rpsi/opt.jData.a;
else
    xvar = opt.jData.rpsi;
end
h = plot(xvar, opt.jData.PI_ASCOT./opt.jData.Qi_QASCOT);
lgd_h = [h];
lgd_txt = {'Experiment (ASCOT)'};

% Plot fluxes from gs2

if ~isempty(opt.gs2_fluxFile)
    if opt.nrm_gs2
        xvar = flx.rpsi/opt.jData.a;
    else
        xvar = flx.rpsi;
    end
    hold on
    lgd_h(end+1) = plot(xvar, flx.PI_gs2.*flx.PINorm./(flx.Qi_gs2.*flx.QNorm), ...
                        'Marker', '.', ...
                        'MarkerSize', 20);
    lgd_txt{end+1} = 'GS2';
end

% Add experimental confidence area

if opt.nrm_gs2
    xvar = opt.jData.rpsi/opt.jData.a;
else
    xvar = opt.jData.rpsi;
end
color = get(h, 'Color');
alpha = 0.3;
hold on
confid_area(gcf, xvar, opt.jData.PI_ASCOT./opt.jData.Qi_QASCOT*0.8, ...
    opt.jData.PI_ASCOT./opt.jData.Qi_QASCOT*1.2, color, alpha)

% Fine-tune figure

if opt.showTitle
    ttl = ['shot=' num2str(opt.jData.shot) ...
        ', ijp=' num2str(opt.jData.ijp) ...
        ', idxAscot=' num2str(opt.jData.idxAscot)];
    title(ttl, 'FontSize',16)
end
grid on
if opt.nrm_gs2
    xlab = '$r_\psi/a$';
else
    xlab = '$r_\psi$ [m]';
end
xlabel(xlab)
ylabel('$\sum_s \Pi_s / Q_i$')
if opt.showAllCodes || ~isempty(opt.gs2_fluxFile)
    legend(lgd_h, lgd_txt, 'Location', 'NorthEast','FontSize',14)
end
if ~isempty(opt.ylim_ratio)
    ylim(opt.ylim_ratio)
elseif min(opt.jData.PI_ASCOT./opt.jData.Qi_QASCOT) > 0
    ylim([0 inf])
end

end
