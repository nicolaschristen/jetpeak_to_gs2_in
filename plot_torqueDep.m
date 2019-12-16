%% Plot experimental torque deposition profiles obtained from ASCOT.
% The torque deposition is expressed in [Pa].
%
% Input :   ijp -- vector of JETPEAK indices
%
% Output:   -
%
function plot_torqueDep(ijp)

for ishot = 1:numel(ijp)

    jData = read_jData(ijp(ishot));

    if ~isempty(jData.idxAscot)
    
        % Total source of torque [Pa]
        S_PEN = jData.srcL_PENCIL;
        S_ASC = jData.srcL_ASCOT;
        
        ttl = ['shot=' num2str(jData.shot) ...
            ', ijp=' num2str(jData.ijp) ...
            ', idxAscot=' num2str(jData.idxAscot)];

        figure

        h = plot(jData.rpsi/jData.a, S_ASC);
        lgd = {'ASCOT'};
        color = get(h,'Color');
        if ~isnan(S_PEN(1))
            hold on
            plot(jData.rpsi/jData.a, S_PEN, 'LineStyle', '--', 'Color', color)
            lgd{end+1} = 'PENCIL';
        end
        grid on
        xlabel('$r_\psi/a$')
        ylabel('$S_{\Pi}$ [Pa]')
        legend(lgd,'Location','SouthEast')
        title(ttl,'FontSize',16)

        fprintf('\nPlotted ijp=%d\n\n', ijp(ishot))

    else

        fprintf('\nJPI=%d has no torque deposition data.\n\n', jData.ijp)

    end

end

end
