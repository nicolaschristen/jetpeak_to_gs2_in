%% Compute GS2 parameters mach and g_exb
%
% Input :   ijp -- shot index in JETPEAK
%           itransp -- shot index in TRANSP
%           psinrm_in -- sqrt(psipol/psipol_LCFS) at flux surf. of interest
%           TRANSP -- structure obtained from loading TRANSP DB
%           EFIT -- structure obtained from loading JETPEAK DB
%           ION -- structure obtained from loading JETPEAK DB
%
% Output:   mach -- Omega_phi * rref/vtref
%           g_exb -- dOmega_phi/dr_psi * r_psi/q * rref/vtref
%
function [mach,g_exb] = compute_mach_g_exb(ijp,itransp,psinrm_in,TRANSP,EFIT,ION)

% Need access to ../interpol.m
addpath(genpath('..'))

% Elementary charge
e=1.602e-19;

% Poloidal flux (GS2 definition, i.e. without toroidal 2pi factor)
psiflu=permute(TRANSP.T.PLFLX,[2 1]); % on flux surface grid, [T m^2]
psiflu=psiflu(:,itransp);
nflxsurf=numel(psiflu); % number of flux surf.
[~,iflxsurf]=min(abs(psiflu-psinrm_in^2*psiflu(end))); % flux-grid point closest to user specified radial location

% Normalized sqrt(psi)
sqrt_psin_TRANSP=zeros(1,nflxsurf);
for i=1:nflxsurf
    sqrt_psin_TRANSP(i)=sqrt(TRANSP.T.PLFLX(itransp,i)/TRANSP.T.PLFLX(itransp,end));
end

% Radial coordinates in meters
a=(TRANSP.G.RMAJM(itransp,end)-TRANSP.G.RMAJM(itransp,1))/2.; % GS2 Lref
rpsi_TRANSP=zeros(1,nflxsurf); % GS2 definition of rho for irho=2, not yet normalized
for indx=1:nflxsurf
    rpsi_TRANSP(indx)=...
        (TRANSP.G.RMAJM(itransp,nflxsurf+1+indx) ...
          - TRANSP.G.RMAJM(itransp,nflxsurf+1-indx))/2.;
end

%% Extract profiles from EFIT, EL and ION

% Normalized sqrt(psi) grid used in chain2
sqrt_psin_chain2=linspace(0.,1.,21);

% Safety factor
q=interpol(EFIT.RMJO(ijp,:),EFIT.Q(ijp,:),TRANSP.G.RMAJM(itransp,nflxsurf+2:end));

% Temperature [eV]
ti=interpol(sqrt_psin_chain2,ION.TI(ijp,:),sqrt_psin_TRANSP);

% Tor. angular vel.
omega=interpol(sqrt_psin_chain2,ION.ANGF(ijp,:),sqrt_psin_TRANSP);
domega_drho=interpol(rpsi_TRANSP,omega,rpsi_TRANSP,1);

% Reference quantities
tref=ti(iflxsurf);
mref=2.*1.67e-27;
vtref=sqrt(2.*e*tref/mref);

% dist_fn_knobs
mach=omega(iflxsurf)*a/vtref;
g_exb=a/vtref*rpsi_TRANSP(iflxsurf)/q(iflxsurf)*domega_drho(iflxsurf);

end