%% Generate NBI deposition profiles of heat and torque
% by modifying manually-fitted experimental profiles.
% A skewed Gaussian function models the shape of profiles.
%
% The user can choose to set:
%
%   -- the beam power with cP = P/P_expt, [0;+inf[
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
%         fit -- structure with fitting parameters
%                (width & nrm) of the manual fit for
%                the experimental profiles of each
%                source term (SQi, SQe and SPi),
%                e.g. fit.SQi.width = 0.35
%         params -- structure with user-set parameters
%                   (cP, cl and cE) for the modified
%                   profiles, e.g. params.cP = 0.5
%
% Output: user -- structure of the user-specified
%                 profiles. Contains mainly srcQi,
%                 srcQe, srcPi, rpsi and a (the
%                 normalising length in GS2).
%         expt -- same as user, but for the manually-
%                 fitted experimental profiles.
%
% NB: The experimental profiles are assumed to peak
% at the tangential flux surface (i.e. cE=0.5), which
% is assumed to be the magnetic axis (i.e. cl=1). They
% are also assumed to have been fitted with a
% skewness of zero.
%
function [user, expt] = set_userDepo(ijp, params, fit, varargin)


% Read optional input arguments
options_default = struct( 'jData', [], ...
                          'plotverbose', 1, ...
                          'nrm_gs2', 1 );
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
expt.SQi.width = fit.SQi.width;
user.SQi.width = expt.SQi.width;

% Electron heating source
expt.SQe.width = fit.SQe.width;
user.SQe.width = expt.SQe.width;

% Total torque source
expt.SPi.width = fit.SPi.width;
user.SPi.width = expt.SPi.width;

%    ------------    %

% Over-all scaling factor of the
% manually-fitted experimental profile

% Ion heating source
expt.SQi.nrm = fit.SQi.nrm;

% Electron heating source
expt.SQe.nrm = fit.SQe.nrm;

% Total torque source
expt.SPi.nrm = fit.SPi.nrm;

%    ------------    %

% Deposition params for the power of the beam

expt.cP = 1;
user.cP = params.cP;

%    ------------    %

% Deposition params for the launching angle lambda

% Assumed experimental value (ie tangent to magnetic axis):
expt.cl = 1;
user.cl = params.cl;
% Determine launching angle at which a beam would
% become tangent to the magnetic axis
lMag = get_lMag(jData);
% Determine the flux surface to which the user-specified
% beam trajectory is tangent
rtan = get_rtan(jData, user.cl);

%    ------------    %

% Deposition params for the energy of the beam

% Maximum possible value:
cEMax = 1;
% Assumed experimental value (ie beam just about reaches the tangent point):
expt.cE = cEMax/2;
user.cE = params.cE;

%    ------------    %

% Set the skewness

expt.SQi.skew = 0;
expt.SQe.skew = 0;
expt.SPi.skew = 0;
S = 10;
if user.cE == cEMax/2
    % Beam just about reaches rtan
    user.SQi.skew = 0;
elseif user.cE < cEMax/2
    % Beam undershoots
    user.SQi.skew = S;
else
    % Beam overshoots
    user.SQi.skew = -S;
end
% Set other quantities to be the same
user.SQe.skew = user.SQi.skew;
user.SPi.skew = user.SQi.skew;

%    ------------    %

% Determine the location of the peak in the deposition

expt.SQi.mpos = 0;
expt.SQe.mpos = 0;
expt.SPi.mpos = 0;
% linear decrease when undershooting,
% and increasing back when overshooting
user.SQi.mpos = (max(jData.rpsi)-rtan)/(cEMax/2) ...
                   * abs(user.cE-cEMax/2) + rtan;
% Set other quantities to be the same
user.SQe.mpos = user.SQi.mpos;
user.SPi.mpos = user.SQi.mpos;

%    ------------    %

% Determine the overall scaling factor to obtain
% a beam with the user-specified power to the ions

% TODO: should we add the electrons when setting the power ???
[expt.P, user.P, user.SQi.nrm] = get_P_and_nrm(jData, expt.SQi, user.SQi, user.cP);

% Then adapt normalisations of all other quantities accordingly
factor = user.SQi.nrm/expt.SQi.nrm;
user.SQe.nrm = factor * expt.SQe.nrm;
user.SPi.nrm = factor * expt.SPi.nrm;

%    ------------    %

% Produce user-specified deposition profiles

user.rpsi = jData.rpsi;
user.a = jData.a;
user.srcQi = gauSkew(jData.rpsi, user.SQi.mpos, user.SQi.width, user.SQi.nrm, user.SQi.skew);
user.srcQe = gauSkew(jData.rpsi, user.SQe.mpos, user.SQe.width, user.SQe.nrm, user.SQe.skew);
user.srcPi = gauSkew(jData.rpsi, user.SPi.mpos, user.SPi.width, user.SPi.nrm, user.SPi.skew);

% And store the manually fitted curves

expt.rpsi = jData.rpsi;
expt.a = jData.a;
expt.srcQi = gauSkew(jData.rpsi, expt.SQi.mpos, expt.SQi.width, expt.SQi.nrm, expt.SQi.skew);
expt.srcQe = gauSkew(jData.rpsi, expt.SQe.mpos, expt.SQe.width, expt.SQe.nrm, expt.SQe.skew);
expt.srcPi = gauSkew(jData.rpsi, expt.SPi.mpos, expt.SPi.width, expt.SPi.nrm, expt.SPi.skew);

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
    lgd_user = [ '$c_P=$' num2str(user.cP) ', $c_\lambda=$' ...
        num2str(user.cl) ', $c_E=$' num2str(user.cE) ];

    figure
    plot(xvar,expt.srcQi)
    hold on
    plot(xvar,user.srcQi)
    xlabel(xlab)
    ylabel('$S_{Q,i}$ [W/m$^3$]')
    legend('Fit to experiment', lgd_user)
    grid on

end

end % end of function


%    ============    %
%    ============    %


function [P_expt, P_user, nrm_user] = get_P_and_nrm(jData, params_expt, params_user, cP)

% First compute the power deposited in the experiment
intgrd = @(r) gauSkew(r, params_expt.mpos, params_expt.width, 1, params_expt.skew);
P_expt = params_expt.nrm * vol_integral(jData, intgrd);

% Then compute the normalisation needed by the user-specified
% profile to get the requested power.
intgrd = @(r) gauSkew(r, params_user.mpos, params_user.width, 1, params_user.skew);
nrm_user = P_expt * cP / vol_integral(jData, intgrd);
P_user = nrm_user * vol_integral(jData, intgrd);

end


%    ============    %
%    ============    %


function intgrl = vol_integral(jData, radial_fun)

intgrl = sum(jData.dV .* radial_fun(jData.rpsi));

end
