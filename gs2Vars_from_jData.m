%% Compute GS2 input parameters from JETPEAK/TRANSP data
%
% Input :   jData -- structure array with dimensionful plasma parameters
%
% Output:   gs2Vars -- structure with radial profiles
%                      of dimensionless GS2 parameters
%
% All quantities are in SI units,
% except for temperatures which are in eV
%
function gs2Vars = gs2Vars_from_jData(jData)

% mag. const. in SI units
mu0=pi*4.e-7;
% elementary charge
e=1.602e-19;

% Reference quantities
a=jData.a;
Bref=jData.Bref;
nref=jData.nref;
tref=jData.tref;
mref=jData.mref;
vthref=jData.vthref;
rhoref=jData.rhoref;
rhostar=jData.rhostar;

% theta_grid_parameters
gs2Vars.Rmaj=jData.Rmaj/a;
gs2Vars.R_geo=jData.Rgeo/a;
gs2Vars.rhoc=jData.rpsi/a;
gs2Vars.shift=jData.dR_drpsi;
gs2Vars.qinp=jData.q;
gs2Vars.shat=jData.shat;
gs2Vars.akappa=jData.kappa;
gs2Vars.akappri=a*jData.dkappa_drpsi;
gs2Vars.tri=jData.delta;
gs2Vars.tripri=a*jData.ddelta_drpsi;

% parameters
gs2Vars.beta=nref*e.*tref./(Bref.^2/(2.*mu0));
gs2Vars.zeff=jData.Zeff;

% theta_grid_eik_knobs
gs2Vars.beta_prime_input=gs2Vars.beta./(nref*e.*tref)*a.*jData.dp_drpsi;

% dist_fn_knobs
gs2Vars.mach=jData.omega*a./vthref;
gs2Vars.g_exb=a./vthref.*jData.rpsi./jData.q.*jData.domega_drpsi;

% Deuterium ions 
% species_parameters_1
gs2Vars.tprim1=-1.*a./jData.ti.*jData.dti_drpsi;
gs2Vars.fprim1=-1.*a./jData.ni.*jData.dni_drpsi;
gs2Vars.dens1=jData.ni./nref;
gs2Vars.temp1=jData.ti./tref;
gs2Vars.vnewk1=jData.nu_ii*a./vthref;

% Electrons
% species_parameters_2
gs2Vars.tprim2=-1.*a./jData.te.*jData.dte_drpsi;
gs2Vars.fprim2=-1.*a./jData.ne.*jData.dne_drpsi;
gs2Vars.dens2=jData.ne./nref;
gs2Vars.temp2=jData.te./tref;
gs2Vars.vnewk2=jData.nu_ee*a./vthref;

% Carbon
% species_parameters_3
gs2Vars.tprim3=-1.*a./jData.tc.*jData.dtc_drpsi;
gs2Vars.fprim3=-1.*a./jData.nc.*jData.dnc_drpsi;
gs2Vars.dens3=jData.nc./nref;
gs2Vars.temp3=jData.tc./tref;
gs2Vars.vnewk3=jData.nu_cc*a./vthref;

end
