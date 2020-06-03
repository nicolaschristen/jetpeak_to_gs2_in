%% Plotting the power and torque deposition profiles
%
% Input :   ijp -- index of the JETPEAK shot
%           jData -- [kw, []] structure of pre-read data from JETPEAK
%           origParams -- structure with fitting parameters
%                  (width & nrm) of the manual fit for
%                  the experimental profiles of each
%                  source term (SQi, SQe and SPi),
%                  e.g. origParams.SQi.width = 0.35
%                  See set_userDepo for details.
%           usrParams -- cell array of structures containing
%                     user-set parameters (cP, cl and cE) for the modified
%                     profiles, e.g. usrParams.cP = 0.5. See set_userDepo for details.
%           usrProfs -- [kw, []] pre-computed cell array with fields srcQi & srcPI
%                       corresponding to elements of usrParams
%           nrm_gs2 -- [kw, 0] normalise plotted quantities to GS2 units
%
% Output:   -
%           
function plot_depo(ijp, varargin)


% Read optional input arguments
options_default = struct( 'jData', [], ...
                          'usrParams', [], ...
                          'nrm_gs2', 0, ...
                          'origParams', [], ...
                          'usr_profs', [] );
opt = get_optargin(options_default, varargin);

% Only read jData if it is not provided by the user
if isempty(opt.jData)
    jData = read_jData(ijp);
else
    jData = opt.jData;
end

% Get number of user-specified profiles
nProf = numel(opt.usrParams);

% Apply nomrmalisation to x-axis
if opt.nrm_gs2
    xvar = jData.rpsi/jData.a;
    xlab = '$r_\psi/a$';
else
    xvar = jData.rpsi;
    xlab = '$r_\psi$ [m]';
end


%    ------------    %

% Plot the power deposition profile

figure

% Experimental
if opt.nrm_gs2
    yvar = []; % TODO
    ylab = ''; % TODO
else
    yvar = jData.srcE_i_QASCOT;
    ylab = '$S_{Q,i}$ [W/m$^3$]';
end
lgd_h = [plot(xvar, yvar)];
lgd_txt = {'Experiment (ASCOT)'};
xlabel(xlab)
ylabel(ylab)
ylim([0, inf])
grid on

% Show skewed Gaussian fit
if ~isempty(opt.origParams)
    if opt.nrm_gs2
        yvar = []; % TODO
    else
        yvar = gauSkew( jData.rpsi, opt.origParams.SQi.mpos, opt.origParams.SQi.width, ...
                        opt.origParams.SQi.nrm, opt.origParams.SQi.skew );
    end
    hold on
    lgd_h(end+1) = plot(xvar, yvar, 'k-', 'LineWidth', 1);
    lgd_txt{end+1} = 'Fit (Gaussian)';
end

% User-specified
for iProf = 1:nProf
    if opt.nrm_gs2
        yvar = []; % TODO
    else
        if isempty(opt.usrProfs)
            [usrProfs, ~] = set_userDepo(ijp, opt.origParams, opt.usrParams{iProf}, ...
                                     'jData', jData);
            yvar = usrProfs.srcQi;
        else
            yvar = opt.usrProfs{iProf}.srcQi;
        end
    end
    hold on
    lgd_h(end+1) = plot(xvar, yvar);
    lgd_txt{end+1} = [ '$c_P=$' num2str(opt.usrParams{iProf}.cP) ', $c_\lambda=$' ...
        num2str(opt.usrParams{iProf}.cl) ', $c_E=$' num2str(opt.usrParams{iProf}.cE) ];
end
    
legend(lgd_h, lgd_txt, 'Location', 'NorthEast','FontSize',14)


%    ------------    %

% Plot the torque deposition profile

figure

% Experimental
if opt.nrm_gs2
    yvar = []; % TODO
    ylab = ''; % TODO
else
    yvar = jData.srcL_ASCOT;
    ylab = '$S_{\Pi}$ [Pa]';
end
lgd_h = [plot(xvar, yvar)];
lgd_txt = {'Experiment (ASCOT)'};
xlabel(xlab)
ylabel(ylab)
ylim([0, inf])
grid on

% Show skewed Gaussian fit
if ~isempty(opt.origParams)
    if opt.nrm_gs2
        yvar = []; % TODO
    else
        yvar = gauSkew( jData.rpsi, opt.origParams.SPi.mpos, opt.origParams.SPi.width, ...
                        opt.origParams.SPi.nrm, opt.origParams.SPi.skew );
    end
    hold on
    lgd_h(end+1) = plot(xvar, yvar, 'k-', 'LineWidth', 1);
    lgd_txt{end+1} = 'Fit (Gaussian)';
end

% User-specified
for iProf = 1:nProf
    if opt.nrm_gs2
        yvar = []; % TODO
    else
        yvar = opt.usrParams{iProf}.srcPI;
        if isempty(opt.usrProfs)
            [usrProfs, ~] = set_userDepo(ijp, opt.origParams, opt.usrParams{iProf}, ...
                                     'jData', jData);
            yvar = usrProfs.srcPI;
        else
            yvar = opt.usrProfs{iProf}.srcPI;
        end
    end
    hold on
    lgd_h(end+1) = plot(xvar, yvar);
    lgd_txt{end+1} = [ '$c_P=$' num2str(opt.usrParams{iProf}.cP) ', $c_\lambda=$' ...
        num2str(opt.usrParams{iProf}.cl) ', $c_E=$' num2str(opt.usrParams{iProf}.cE) ];
end
    
legend(lgd_h, lgd_txt, 'Location', 'NorthEast','FontSize',14)

end
