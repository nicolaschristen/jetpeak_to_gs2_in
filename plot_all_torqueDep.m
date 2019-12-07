%% Plot experimental torque deposition profiles obtained from ASCOT.
% The torque deposition is expressed in [Pa].
%
% Input :   ijp -- vector of JETPEAK indices of shot
%           sampAscot -- vector of corresponding ASCOT samples
%
% Output:   -
%
function plot_all_torqueDep(ijp,sampAscot)

load('~/codes/jetpeak_and_gs2/databases/ASCOT_2019_11.mat','ASCOT')

figure
labelTxt = {};

for ishot = 1:numel(ijp)

    idxAscot = find(ASCOT.RUNDIR==sampAscot(ishot));

    skip_I = 1;
    jData = read_jData(ijp(ishot), skip_I);
    
    % Total source of torque [Pa]
    S_Pi = ASCOT.COLLTQI(idxAscot,:) + ASCOT.COLLTQE(idxAscot,:) ...
        + ASCOT.COLLTQIMP(idxAscot,:) + ASCOT.JXBTORQ(idxAscot,:) ...
        - ASCOT.TQTH(idxAscot,:);

    % Interpolate from CHAIN2 to TRANSP radial grid
    S_Pi = interpol(jData.sqrt_psin_chain2,S_Pi,jData.sqrt_psin_TRANSP);
    
    labelTxt{end+1} = ['shot=' num2str(jData.shot) ', ' ...
        'ijp=' num2str(ijp(ishot)) ', ' ...
        'idxAscot=' num2str(sampAscot(ishot))];
    plot(jData.rpsi/jData.a, S_Pi);
    hold on

end % loop over shots

grid on
xlabel('$r_\psi/a$')
ylabel('$S_{\Pi}$ [Pa]')
%legend(lgd_h, lgd_txt, 'Location', 'NorthEast','FontSize',14)

% The user can highlight a line by left-clicking on it.
user_select_line(labelTxt)

end
