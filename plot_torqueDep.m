%% Plot experimental torque deposition profiles obtained from ASCOT.
% The torque deposition is expressed in [Pa].
%
% Input :   ijp -- vector of JETPEAK indices
%           nrm_gs2 -- [optional, =0] use GS2 normalisation
%                      for x-axis
%
% Output:   -
%
function plot_torqueDep(ijp, varargin)


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
    
        % Total source of torque [Pa]
        S_PEN = jData.srcL_PENCIL;
        S_ASC = jData.srcL_ASCOT;
        
        ttl = ['shot=' num2str(jData.shot) ...
            ', ijp=' num2str(jData.ijp) ...
            ', idxAscot=' num2str(jData.idxAscot)];

        figure

        h = plot(xvar, S_ASC);
        lgd = {'ASCOT'};
        color = get(h,'Color');
        if ~isnan(S_PEN(1))
            hold on
            plot(xvar, S_PEN, 'LineStyle', '--', 'Color', color)
            lgd{end+1} = 'PENCIL';
        end
        grid on
        xlabel(xlab)
        ylabel('$S_{\Pi}$ [Pa]')
        legend(lgd,'Location','SouthEast')
        title(ttl,'FontSize',16)

        fprintf('\nPlotted ijp=%d\n\n', ijp(ishot))

    else

        fprintf('\nJPI=%d has no torque deposition data.\n\n', jData.ijp)

    end

end

end
