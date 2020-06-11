%% Plot experimental momentum flux obtained from ASCOT.
% The fluxes are expressed in GS2 units [nref*mref*vthref^2*Lref*rhostar^2].
%
% Input :   ijp -- JETPEAK index of shot
%           ylim -- [kw, []] ylim for plot
%           jData -- [kw, []] provide pre-read JETPEAK data
%           gs2_fluxFile -- [kw, ''] file containing GS2 fluxes
%                           to be plotted
%           nrm_gs2 -- [kw, 0] normalise plotted quantities to GS2 units
%           showTitle -- [kw, 1] add title to plots
%           showAllCodes -- [kw, 0] plot other deposition codes than
%                           ASCOT, eg PENCIL
%           trinity_norm -- [kw, 0] if true, gs2 flux dotted with gradPsi
%                           else dotted with grad(x).
%
% Output:   r -- dimensionful radial coordinate
%           PI_ASC -- momentum flux from ASCOT
%
function [r, PI_ASC] = plot_momFlux(ijp, varargin)

opt_defaults = struct( 'ylim', [], ...
                       'jData', [], ...
                       'gs2_fluxFile', '', ...
                       'nrm_gs2', 0, ...
                       'showTitle', 1, ...
                       'showAllCodes', 0, ...
                       'trinity_norm', 0 );
opt = get_optargin(opt_defaults, varargin);

% Elementary charge
e=1.602e-19;
    
% Only read data if it has not been passed as an argument
if isempty(opt.jData)
    jData = read_jData(ijp, 'trinity_norm', opt.trinity_norm);
else
    jData = opt.jData;
end

% Output
r = jData.rpsi;
PI_ASC = jData.PI_ASCOT;
    
% Read GS2 fluxes form file

if ~isempty(opt.gs2_fluxFile)
    flx = read_gs2Fluxes( ijp, opt.gs2_fluxFile, 'jData', jData, ...
                          'trinity_norm', opt.trinity_norm );
end

%% Total momentum flux [W]

figure
lgd_h = [];
lgd_txt = {};

% Plot

if opt.nrm_gs2
    xvar = jData.rpsi/jData.a;
    xlab = '$r_\psi/a$';
    ylab = '$\sum_s\Pi_s$ [$n_r$ $r_r$ $m_r$ $v_{thr}^2$ $\rho_\star^2$]';
    yvar_PEN = jData.PI_PENCIL./jData.PINorm;
    yvar_ASC = jData.PI_ASCOT./jData.PINorm;
else
    xvar = jData.rpsi;
    xlab = '$r_\psi$ [m]';
    ylab = '$\sum_s\Pi_s$ [kg/s$^2$]';
    yvar_PEN = jData.PI_PENCIL;
    yvar_ASC = jData.PI_ASCOT;
end

if ~isnan(jData.PI_PENCIL(1)) && opt.showAllCodes

    lgd_h(end+1) = semilogy(xvar, yvar_PEN);
    hold on
    lgd_txt{end+1} = 'Experiment (PENCIL)';
    color = get(lgd_h(end), 'Color');
    alpha = 0.3;
    confid_area(gcf, xvar, yvar_PEN*0.8, ...
        yvar_PEN*1.2, color, alpha)

end

h = semilogy(xvar, yvar_ASC);
lgd_h(end+1) = h;
lgd_txt{end+1} = 'Experiment (ASCOT)';

% Plot flux from  gs2 simulations

if ~isempty(opt.gs2_fluxFile)
    if opt.nrm_gs2
        xvar = flx.rpsi/jData.a;
        yvar = flx.PI_gs2;
    else
        xvar = flx.rpsi;
        yvar = flx.PI_gs2.*flx.PINorm;
    end
    hold on
    lgd_h(end+1) = semilogy(xvar, yvar, ...
                            'Marker', '.', ...
                            'MarkerSize', 20);
    lgd_txt{end+1} = 'GS2';
end

% Add experimental confidence area

if opt.nrm_gs2
    xvar = jData.rpsi/jData.a;
else
    xvar = jData.rpsi;
end
hold on
color = get(h, 'Color');
alpha = 0.3;
confid_area(gcf, xvar, yvar_ASC*0.8, ...
    yvar_ASC*1.2, color, alpha)

% Fine-tune figure

if opt.showTitle
    ttl = ['shot=' num2str(jData.shot) ...
        ', ijp=' num2str(jData.ijp) ...
        ', idxAscot=' num2str(jData.idxAscot)];
    title(ttl, 'FontSize',16)
end

grid on
xlabel(xlab)
ylabel(ylab)
if opt.showAllCodes || ~isempty(opt.gs2_fluxFile)
    legend(lgd_h, lgd_txt, 'Location', 'SouthEast','FontSize',14)
end
if ~isempty(opt.ylim)
    ylim(opt.ylim)
end

end
