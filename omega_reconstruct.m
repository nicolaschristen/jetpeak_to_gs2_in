function omega_reconstruct(ijp, fname, varargin)


% Read optional input arguments
options_default = struct( 'use_Pi_over_Q', 0, ...
                          'depoParams', [] );
opt = get_optargin(options_default, varargin);

%    ------------    %

% Read data for this shot from JETPEAK
jData = read_jData(ijp);

% Either read the fluxes from simulations,
% use them to compute the transport coefficients,
% and save those coefficients to a file.
if isempty(opt.depoParams)
    [flx, tCoeff] = get_transpCoeff_from_gs2Fluxes(ijp, fname, 'jData', jData);
% Or compute the fluxes based on user-specified profiles,
% and read pre-computed transport coefficients from a file.
else
    % TODO
    flx = 0;
    tCoeff = readtable('transpCoeff.xls');
end

%    ------------    %

% Define an rpsi interval on which to reconstruct omega
rpsiIn = min(tCoeff.rpsi);
rpsiOut = max(tCoeff.rpsi);
% Index in jData corresponding to rpsiOut
iout_jData = find(abs(rpsiOut-jData.rpsi) == min(abs(rpsiOut-jData.rpsi)));
nr = 1000;
recon.rpsi = linspace(rpsiIn,rpsiOut,nr);

%    ------------    %

% For quantities determined from simulations, use linear interpolation.
% For user-set or experimental quantities, use splines.
recon.momPinch = interp1(tCoeff.rpsi,tCoeff.momPinch,recon.rpsi);
recon.momDif = interp1(tCoeff.rpsi,tCoeff.momDif,recon.rpsi);
if isempty(opt.depoParams)
    recon.PI = interp1(flx.rpsi,flx.PI_gs2.*flx.PINorm,recon.rpsi);
else
    recon.PI = interpol(flx.rpsi,flx.PI_gs2.*flx.PINorm,recon.rpsi);
end
recon.nref = interpol(jData.rpsi,jData.nref,recon.rpsi);
recon.Rmaj = interpol(jData.rpsi,jData.Rmaj,recon.rpsi);

%    ------------    %

% Integrate inwards using trapezoidal rule

dr = recon.rpsi(2) - recon.rpsi(1);

% First compute exponent
g = zeros(1,nr);
if ~opt.use_Pi_over_Q
    intgrd = recon.momPinch./(recon.momDif.*recon.Rmaj);
else
    % TODO
    intgrd = 0;
end
for ir = nr-1:-1:1
    g(ir) = g(ir+1) + (intgrd(ir)+intgrd(ir+1)) * dr/2.0;
end
% Then compute main integral
f = zeros(1,nr);
if ~opt.use_Pi_over_Q
    intgrd = recon.PI./(recon.momDif.*recon.nref.*jData.mref.*recon.Rmaj.^2) .* exp(-1*g);
else
    % TODO
    intgrd = 0;
end
for ir = nr-1:-1:1
    f(ir) = f(ir+1) + (intgrd(ir)+intgrd(ir+1)) * dr/2.0;
end

%    ------------    %

% Get and save expression for omega(rpsi)

% Edge BC for the angular frequency
omegaOut = abs(jData.omega(iout_jData));

recon.omega = omegaOut * ones(1,nr) + exp(g).*f;

% Get nearest neighbour on recon grid to rpsi's from the simulations
nr_tCoeff = numel(tCoeff.rpsi);
omega = zeros(nr_tCoeff,1);
for ir = 1:nr_tCoeff
    r = tCoeff.rpsi(ir);
    ir_recon2tCoeff = find(abs(recon.rpsi-r) == min(abs(recon.rpsi-r)));
    omega(ir) = recon.omega(ir_recon2tCoeff);
end

% Write transport coefficients to file
% in same folder as file containing fluxes.
rpsi = tCoeff.rpsi;
omRecon = table(rpsi,omega);
[fpath,~,~] = fileparts(fname);
writetable(omRecon, [fpath '/omegaReconstruct.csv']);

%    ------------    %

% Plot omega

figure
h = [plot(omRecon.rpsi,omRecon.omega,'.-','MarkerSize',20,'LineWidth',1)];
lgd = {'reconstructed'};
if isempty(opt.depoParams)
    hold on
    h(end+1) = plot(jData.rpsi,abs(jData.omega));
    lgd{end+1} = 'experiment';
end
legend(h,lgd)
xlabel('$r_\psi / a$')
ylabel('$\Omega_\phi$ [rad/s]')
grid on

end
