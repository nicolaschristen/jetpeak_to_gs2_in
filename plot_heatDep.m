%% Plot experimental heat deposition profiles obtained from ASCOT.
% The heat deposition is expressed in [W/m^3].
%
% Input :   ijp -- vector of JETPEAK indices
%           nrm_gs2 -- [optional, =0] use GS2 normalisation
%                      for x-axis
%
% Output:   -
%
function plot_heatDep(ijp, varargin)


% Read optional input arguments
options_default = struct( 'nrm_gs2', 0 );
opt = get_optargin(options_default, varargin);

nrm_gs2 = opt.nrm_gs2;

for ishot = 1:numel(ijp)

    jData = read_jData(ijp(ishot));

    if nrm_gs2
        xvar = jData.rpsi/jData.a;
        xlab = '$r_\psi/a$';
    else
        xvar = jData.rpsi;
        xlab = '$r_\psi$ [m]';
    end

    if ~isempty(jData.idxAscot)
    
        % Total source of heat [W/m^3]
        Si_PEN = jData.srcE_i_PENCIL;
        Se_PEN = jData.srcE_e_PENCIL;
        Si_ASC = jData.srcE_i_QASCOT;
        Se_ASC = jData.srcE_e_QASCOT;
        
        ttl = ['shot=' num2str(jData.shot) ...
            ', ijp=' num2str(jData.ijp) ...
            ', idxAscot=' num2str(jData.idxAscot)];

        % Start by plotting ions

        figure

        h = plot(xvar, Si_ASC);
        lgd = {'ASCOT'};
        color = get(h,'Color');
        if ~isnan(Si_PEN(1))
            hold on
            plot(xvar, Si_PEN, 'LineStyle', '--', 'Color', color)
            lgd{end+1} = 'PENCIL';
        end
        grid on
        xlabel(xlab)
        ylabel('$S_{Q,i}$ [W/$m^3$]')
        legend(lgd,'Location','SouthEast')
        title(ttl,'FontSize',16)

        % Then plot the electrons

        figure

        h = plot(xvar, Se_ASC);
        lgd = {'ASCOT'};
        color = get(h,'Color');
        if ~isnan(Se_PEN(1))
            hold on
            plot(xvar, Se_PEN, 'LineStyle', '--', 'Color', color)
            lgd{end+1} = 'PENCIL';
        end
        grid on
        xlabel(xlab)
        ylabel('$S_{Q,e}$ [W/$m^3$]')
        legend(lgd,'Location','SouthEast')
        title(ttl,'FontSize',16)


        fprintf('\nPlotted ijp=%d\n\n', ijp(ishot))

    else

        fprintf('\nJPI=%d has no torque deposition data.\n\n', jData.ijp)

    end

end

end
