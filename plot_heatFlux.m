%% Plot experimental heat flux obtained from power balance.
% The fluxes are expressed in GS2 units [nref*Tref*vthref*rhostar^2].
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
%           Qi_ASC -- ion heat flux from ASCOT
%           Qe_ASC -- electron heat flux from ASCOT
%
function [r, Qi_ASC, Qe_ASC] = plot_heatFlux(ijp, varargin)

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

% Function output
r = jData.rpsi;
Qi_ASC = jData.Qi_QASCOT;
Qe_ASC = jData.Qe_QASCOT;
    
% Read GS2 fluxes form file

if ~isempty(opt.gs2_fluxFile)
    flx = read_gs2Fluxes(ijp, opt.gs2_fluxFile, 'jData', jData, ...
                         'trinity_norm', opt.trinity_norm );
end

%% Net ion heat flux [W]

figure
lgd_h = [];
lgd_txt = {};

% Plot

if opt.nrm_gs2
    xvar = jData.rpsi/jData.a;
    xlab = '$r_\psi/a$';
    ylab = '$Q_i$ [$n_r T_r v_{thr} \rho_\star^2$]';
    yvar_PEN = jData.Qi_PENCIL./jData.QNorm;
    yvar_ASC = jData.Qi_QASCOT./jData.QNorm;
else
    xvar = jData.rpsi;
    xlab = '$r_\psi$ [m]';
    ylab = '$Q_i$ [W/m$^2$]';
    yvar_PEN = jData.Qi_PENCIL;
    yvar_ASC = jData.Qi_QASCOT;
end

if ~isnan(jData.Qi_PENCIL(1)) && opt.showAllCodes

    h = semilogy(xvar, yvar_PEN);
    lgd_h(end+1) = h;
    hold on
    lgd_txt{end+1} = 'PENCIL';
    color = get(h, 'Color');
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
        yvar = flx.Qi_gs2;
    else
        xvar = flx.rpsi;
        yvar = flx.Qi_gs2.*flx.QNorm;
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



%% Net electron heat flux [W]

figure
lgd_h = [];
lgd_txt = {};

% Plot

if opt.nrm_gs2
    xvar = jData.rpsi/jData.a;
    xlab = '$r_\psi/a$';
    ylab = '$Q_e$ [$n_r T_r v_{thr} \rho_\star^2$]';
    yvar_PEN = jData.Qe_PENCIL./jData.QNorm;
    yvar_ASC = jData.Qe_QASCOT./jData.QNorm;
else
    xvar = jData.rpsi;
    xlab = '$r_\psi$ [m]';
    ylab = '$Q_e$ [W/m$^2$]';
    yvar_PEN = jData.Qe_PENCIL;
    yvar_ASC = jData.Qe_QASCOT;
end

if ~isnan(jData.Qe_PENCIL(1)) && opt.showAllCodes

    h = semilogy(xvar, yvar_PEN);
    lgd_h(end+1) = h;
    hold on
    lgd_txt{end+1} = 'PENCIL';
    color = get(h, 'Color');
    alpha = 0.3;
    confid_area(gcf, xvar, yvar_PEN*0.8, ...
        yvar_PEN*1.2, color, alpha)

end

h = semilogy(xvar, yvar_ASC);
lgd_h(end+1) = h;
lgd_txt{end+1} = 'Experiment (ASCOT)';

% Plot flux from gs2 simulations

if ~isempty(opt.gs2_fluxFile) && iscol(flx,'Qe_gs2')
    if opt.nrm_gs2
        xvar = flx.rpsi/jData.a;
        yvar = flx.Qe_gs2;
    else
        xvar = flx.rpsi;
        yvar = flx.Qe_gs2.*flx.QNorm;
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
if opt.showAllCodes || (~isempty(opt.gs2_fluxFile) && isfield(flx,'Qe_gs2'))
    legend(lgd_h, lgd_txt, 'Location', 'SouthEast','FontSize',14)
end
if ~isempty(opt.ylim)
    ylim(opt.ylim)
end

end
