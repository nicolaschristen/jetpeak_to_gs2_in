%% Plot experimental heat deposition profiles obtained from ASCOT.
% The heat deposition is expressed in [W/m^3].
%
% Input :   ijp -- vector of JETPEAK indices
%
% Output:   -
%
function plot_heatDep(ijp)

for ishot = 1:numel(ijp)

    jData = read_jData(ijp(ishot));

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

        h = plot(jData.rpsi/jData.a, Si_ASC);
        lgd = {'ASCOT'};
        color = get(h,'Color');
        if ~isnan(Si_PEN(1))
            hold on
            plot(jData.rpsi/jData.a, Si_PEN, 'LineStyle', '--', 'Color', color)
            lgd{end+1} = 'PENCIL';
        end
        grid on
        xlabel('$r_\psi/a$')
        ylabel('$S_{Q,i}$ [W/$m^3$]')
        legend(lgd,'Location','SouthEast')
        title(ttl,'FontSize',16)

        % Then plot the electrons

        figure

        h = plot(jData.rpsi/jData.a, Se_ASC);
        lgd = {'ASCOT'};
        color = get(h,'Color');
        if ~isnan(Se_PEN(1))
            hold on
            plot(jData.rpsi/jData.a, Se_PEN, 'LineStyle', '--', 'Color', color)
            lgd{end+1} = 'PENCIL';
        end
        grid on
        xlabel('$r_\psi/a$')
        ylabel('$S_{Q,e}$ [W/$m^3$]')
        legend(lgd,'Location','SouthEast')
        title(ttl,'FontSize',16)


        fprintf('\nPlotted ijp=%d\n\n', ijp(ishot))

    else

        fprintf('\nJPI=%d has no torque deposition data.\n\n', jData.ijp)

    end

end

end
