function fluxes = get_fluxes_from_userDepo(params, fit)

% Keep same width for the user-specified profile
% as in the manually-fitted experimental one.
expt.width = fit.width;
user.width = expt.width;

%    ------------    %

% Read over-all scaling factor of the
% manually-fitted experimental profile
expt.nrm = fit.nrm;

%    ------------    %

% Deposition params for the power of the beam

% Experimental value:
expt.cP = 1;
user.cP = params.cP;

%    ------------    %

% Deposition params for the launching angle lambda

% Assumed experimental value:
expt.cl = 1;
user.cl = params.cl;
% Determine launching angle at which a beam would
% become tangent to the magnetic axis
lMag = TODO;
% Determine the flux surface to which the
% beam trajectory is tangent
rtan = TODO;

%    ------------    %

% Deposition params for the energy of the beam

% Maximum possible value:
cEMax = 1;
% Assumed experimental value:
expt.cE = cEMax/2;
user.cE = params.cE;

%    ------------    %

% Set the skewness

expt.skew = 0;
S = 10;
if user.cE = cEMax/2
    % Beam just about reaches rtan
    user.skew = 0;
else if user.cE < cEMax/2
    % Beam undershoots
    user.skew = S;
else
    % Beam overshoots
    user.skew = -S;
end

%    ------------    %

% Determine the location of the peak in the deposition

expt.mpos = 0;
user.mpos = (rmax-rtan)/(cEMax/2) * abs(user.cE-cEMax/2) + rtan;

%    ------------    %

% Determine the overall scaling factor to obtain
% a beam with the user-specified power

expt.P = TODO;
user.P = user.cP*expt.P;
user.int_noNorm = TODO;
user.nrm = user.P/user.int_noNorm;

%    ------------    %

% Produce the user-specified deposition profiles

TODO

%    ------------    %

% Compute the fluxes associated with these deposition profiles

TODO


end
