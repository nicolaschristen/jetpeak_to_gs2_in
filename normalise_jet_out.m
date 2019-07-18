%% Apply GS2 normalisations to parameters extracted from JETPEAK/TRANSP
%
% Input :   jet_out -- structure array with dimensionful plasma parameters
%           add_carbon -- =0 neglects C and adapts n_ion so that QN holds
%                         =1 adds Carbon as an impurity
%
% Output:   gs2_in -- structure with dimensionless GS2 parameters
%
% All quantities are in SI units,
% except for temperatures which are in eV
%
function gs2_in=normalise_jet_out(jet_out, add_carbon)

% mag. const. in SI units
mu0=pi*4.e-7;
% elementary charge
e=1.602e-19;

% Reference quantities
a=jet_out.a;
Bref=jet_out.Bref;
nref=jet_out.nref;
tref=jet_out.tref;
mref=jet_out.mref;
vtref=sqrt(2.*e*tref/mref);

% theta_grid_parameters
gs2_in.Rmaj=jet_out.Rmaj/a;
gs2_in.R_geo=jet_out.Rgeo/a;
gs2_in.rhoc=jet_out.rho/a;
gs2_in.shift=jet_out.dR_drho;
gs2_in.qinp=jet_out.q;
gs2_in.shat=jet_out.rho/jet_out.q * jet_out.dq_drho;
gs2_in.akappa=jet_out.kappa;
gs2_in.akappri=a*jet_out.dkappa_drho;
gs2_in.tri=jet_out.delta;
gs2_in.tripri=a*jet_out.ddelta_drho;

% parameters
gs2_in.beta=nref*e*tref/(Bref^2/(2.*mu0));
gs2_in.zeff=jet_out.Zeff;

% theta_grid_eik_knobs
gs2_in.beta_prime_input=gs2_in.beta/jet_out.p*a*jet_out.dp_drho;

% dist_fn_knobs
gs2_in.mach=jet_out.omega*a/vtref;
gs2_in.g_exb=a/vtref*jet_out.rho/jet_out.q*jet_out.domega_drho;

% Deuterium ions 
% species_parameters_1
gs2_in.tprim1=-1.*a/jet_out.ti*jet_out.dti_drho;
gs2_in.fprim1=-1.*a/jet_out.ni*jet_out.dni_drho;
gs2_in.dens1=jet_out.ni/nref;
gs2_in.temp1=jet_out.ti/tref;

% Electrons
% species_parameters_2
gs2_in.tprim2=-1.*a/jet_out.te*jet_out.dte_drho;
gs2_in.fprim2=-1.*a/jet_out.ne*jet_out.dne_drho;
gs2_in.dens2=jet_out.ne/nref;
gs2_in.temp2=jet_out.te/tref;

% Carbon
% species_parameters_3
if add_carbon
    gs2_in.tprim3=-1.*a/jet_out.tc*jet_out.dtc_drho;
    gs2_in.fprim3=-1.*a/jet_out.nc*jet_out.dnc_drho;
    gs2_in.dens3=jet_out.nc/nref;
    gs2_in.temp3=jet_out.tc/tref;
end

end