%% Generate NBI deposition profiles of heat and torque
% by modifying manually-fitted experimental profiles.
% A skewed Gaussian function models the shape of profiles.
%
% The user can choose to set:
%
%   -- the beam power with cP = P/P_origDepo, [0;+inf[
%
%   -- the launching angle with cl = angle/angle_mag,
%      where angle_mag is the angle to become tangent
%      to the magnetic axis, [0;1]
%
%   -- the energy, or penetration of the beam with
%      cE, [0;1], with cE = 0 giving a profile that
%      does not penetrate, 0.5 giving a beam that
%      peaks at the tangent flux surface, and 1.0
%      peaking after traversing the whole plasma.
%
% Input : ijp -- JETPEAK index of the shot
%         origParams -- structure with fitting parameters
%                (width & nrm) of the manual fit for
%                the experimental profiles of each
%                source term (SQi, SQe and SPi),
%                e.g. origParams.SQi.width = 0.35
%         usrParams -- structure with user-set parameters
%                   (cP, cl and cE) for the modified
%                   profiles, e.g. usrParams.cP = 0.5
%
% Output: usrDepo -- structure of the user-specified
%                 profiles. Contains mainly srcQi,
%                 srcQe, srcPi, rpsi and a (the
%                 normalising length in GS2).
%         origDepo -- same as usrDepo, but for the manually-
%                 fitted experimental profiles.
%
% NB: The experimental profiles are assumed to peak
% at the tangential flux surface (i.e. cE=0.5), which
% is assumed to be the magnetic axis (i.e. cl=1). They
% are also assumed to have been fitted with a
% skewness of zero.
%
function [usrDepo, origDepo] = set_userDepo(ijp, origParams, usrParams, varargin)


% Read optional input arguments
options_default = struct( 'jData', [], ...
                          'plotverbose', 0, ...
                          'nrm_gs2', 0 );
opt = get_optargin(options_default, varargin);

% Read JETPEAK data if it has not already been done
if isempty(opt.jData)
    jData = read_jData(ijp);
else
    jData = opt.jData;
end

%    ------------    %

% Keep same width for the user-specified profile
% as in the manually-fitted experimental one.

% Ion heating source
origDepo.SQi.width = origParams.SQi.width;
usrDepo.SQi.width = origDepo.SQi.width;

% Electron heating source
origDepo.SQe.width = origParams.SQe.width;
usrDepo.SQe.width = origDepo.SQe.width;

% Total torque source
origDepo.SPi.width = origParams.SPi.width;
usrDepo.SPi.width = origDepo.SPi.width;

%    ------------    %

% Over-all scaling factor of the
% manually-fitted experimental profile

% Ion heating source
origDepo.SQi.nrm = origParams.SQi.nrm;

% Electron heating source
origDepo.SQe.nrm = origParams.SQe.nrm;

% Total torque source
origDepo.SPi.nrm = origParams.SPi.nrm;

%    ------------    %

% Deposition params for the power of the beam

origDepo.cP = 1;
usrDepo.cP = usrParams.cP;

%    ------------    %

% Deposition params for the launching angle lambda

% Assumed experimental value (ie tangent to magnetic axis):
origDepo.cl = 1;
usrDepo.cl = usrParams.cl;
% Determine launching angle at which a beam would
% become tangent to the magnetic axis
lMag = get_lMag(jData);
% Determine the flux surface to which the user-specified
% beam trajectory is tangent
rtan = get_rtan(jData, usrDepo.cl);

%    ------------    %

% Deposition params for the energy of the beam

% Maximum possible value:
cEMax = 1;
% Assumed experimental value (ie beam just about reaches the tangent point):
origDepo.cE = cEMax/2;
usrDepo.cE = usrParams.cE;

%    ------------    %

% Set the skewness

origDepo.SQi.skew = 0;
origDepo.SQe.skew = 0;
origDepo.SPi.skew = 0;
S = 10;
if usrDepo.cE == cEMax/2
    % Beam just about reaches rtan
    usrDepo.SQi.skew = 0;
elseif usrDepo.cE < cEMax/2
    % Beam undershoots
    usrDepo.SQi.skew = S;
else
    % Beam overshoots
    usrDepo.SQi.skew = -S;
end
% Set other quantities to be the same
usrDepo.SQe.skew = usrDepo.SQi.skew;
usrDepo.SPi.skew = usrDepo.SQi.skew;

%    ------------    %

% Determine the location of the peak in the deposition

origDepo.SQi.mpos = 0;
origDepo.SQe.mpos = 0;
origDepo.SPi.mpos = 0;
% linear decrease when undershooting,
% and increasing back when overshooting
usrDepo.SQi.mpos = (max(jData.rpsi)-rtan)/(cEMax/2) ...
                   * abs(usrDepo.cE-cEMax/2) + rtan;
% Set other quantities to be the same
usrDepo.SQe.mpos = usrDepo.SQi.mpos;
usrDepo.SPi.mpos = usrDepo.SQi.mpos;

%    ------------    %

% Determine the overall scaling factor to obtain
% a beam with the user-specified power to the ions

% TODO: should we add the electrons when setting the power ???
[origDepo.P, usrDepo.P, usrDepo.SQi.nrm] = get_P_and_nrm(jData, origDepo.SQi, usrDepo.SQi, usrDepo.cP);

% Then adapt normalisations of all other quantities accordingly
factor = usrDepo.SQi.nrm/origDepo.SQi.nrm;
usrDepo.SQe.nrm = factor * origDepo.SQe.nrm;
usrDepo.SPi.nrm = factor * origDepo.SPi.nrm;

%    ------------    %

% Produce user-specified deposition profiles

usrDepo.rpsi = jData.rpsi;
usrDepo.a = jData.a;
usrDepo.srcQi = gauSkew(jData.rpsi, usrDepo.SQi.mpos, usrDepo.SQi.width, ...
                        usrDepo.SQi.nrm, usrDepo.SQi.skew);
usrDepo.srcQe = gauSkew(jData.rpsi, usrDepo.SQe.mpos, usrDepo.SQe.width, ...
                        usrDepo.SQe.nrm, usrDepo.SQe.skew);
usrDepo.srcPi = gauSkew(jData.rpsi, usrDepo.SPi.mpos, usrDepo.SPi.width, ...
                        usrDepo.SPi.nrm, usrDepo.SPi.skew);

% And store the manually fitted curves

origDepo.rpsi = jData.rpsi;
origDepo.a = jData.a;
origDepo.srcQi = gauSkew(jData.rpsi, origDepo.SQi.mpos, origDepo.SQi.width, ...
                         origDepo.SQi.nrm, origDepo.SQi.skew);
origDepo.srcQe = gauSkew(jData.rpsi, origDepo.SQe.mpos, origDepo.SQe.width, ...
                         origDepo.SQe.nrm, origDepo.SQe.skew);
origDepo.srcPi = gauSkew(jData.rpsi, origDepo.SPi.mpos, origDepo.SPi.width, ...
                         origDepo.SPi.nrm, origDepo.SPi.skew);

%    ------------    %

% Plotting

if opt.plotverbose

    if opt.nrm_gs2
        xvar = jData.rpsi/jData.a;
        xlab = '$r_\psi/a$';
    else
        xvar = jData.rpsi;
        xlab = '$r_\psi$';
    end

    % Legend for user-specified profiles
    lgd_user = [ '$c_P=$' num2str(usrDepo.cP) ', $c_\lambda=$' ...
        num2str(usrDepo.cl) ', $c_E=$' num2str(usrDepo.cE) ];

    figure
    plot(xvar,origDepo.srcQi)
    hold on
    plot(xvar,usrDepo.srcQi)
    xlabel(xlab)
    ylabel('$S_{Q,i}$ [W/m$^3$]')
    legend('Fit to experiment', lgd_user)
    grid on

end

end % end of function


%    ============    %
%    ============    %


function [P_orig, P_user, nrm_user] = get_P_and_nrm(jData, params_orig, params_user, cP)

% First compute the power deposited in the experiment
intgrd = @(r) gauSkew(r, params_orig.mpos, params_orig.width, 1, params_orig.skew);
P_orig = params_orig.nrm * vol_integral(jData, intgrd);

% Then compute the normalisation needed by the user-specified
% profile to get the requested power.
intgrd = @(r) gauSkew(r, params_user.mpos, params_user.width, 1, params_user.skew);
nrm_user = P_orig * cP / vol_integral(jData, intgrd);
P_user = nrm_user * vol_integral(jData, intgrd);

end


%    ============    %
%    ============    %


function intgrl = vol_integral(jData, radial_fun)

intgrl = sum(jData.dV .* radial_fun(jData.rpsi));

end
