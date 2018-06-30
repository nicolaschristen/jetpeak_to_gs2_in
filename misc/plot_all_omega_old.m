% Elementary charge
e=1.602e-19;

% mag. const. in SI units
mu0=pi*4.e-7;

% Ignore NaN warnings
warning('off','MATLAB:chckxy:IgnoreNaN');

fig_mach_edge=figure;
xlabel('ijp')
ylabel('mach [vthref/a]')
title('GS2 $\omega$ at $\sqrt{\psi_N}=0.8$','Interprester','Latex')

fig_g_exb_edge=figure;
xlabel('ijp')
ylabel('g_exb [q*vthref/(a*rho)]')
title('GS2 $\frac{d\omega}{d\rho}$ at $\sqrt{\psi_N}=0.8$','Interprester','Latex')

fig_mach_in=figure;
xlabel('ijp')
ylabel('mach [vthref/a]')
title('GS2 $\omega$ at $\sqrt{\psi_N}=0.3$','Interprester','Latex')

fig_g_exb_in=figure;
xlabel('ijp')
ylabel('g_exb [q*vthref/(a*rho)]')
title('GS2 $\frac{d\omega}{d\rho}$ at $\sqrt{\psi_N}=0.3$','Interprester','Latex')

for ijp=1:1847
    itransp=find(TRANSP.JPI==ijp);
    if ~isempty(itransp)
        itransp=itransp(1);
        
        % compute GS2 mach & g_exb at the edge
        psinrm_edge=0.8;
        [mach_edge,g_exb_edge]=compute_mach_g_exb(ijp,itransp,psinrm_edge);
        % add point to edge plot
        figure(fig_mach_edge);
        plot(ijp,mach_edge,'bx')
        figure(fig_g_exb_edge);
        plot(ijp,g_exb_edge,'bx')
        
        % same for inner
        psinrm_in=0.3;
        [mach_in,g_exb_in]=compute_mach_g_exb(ijp,itransp,psinrm_in);
        figure(fig_mach_in);
        plot(ijp,mach_in,'bx')
        figure(fig_g_exb_in);
        plot(ijp,g_exb_in,'bx')
              
    end    
end



function [mach,g_exb] = compute_mach_g_exb(ijp,itransp,psinrm_in)

% Poloidal flux (GS2 definition, i.e. without toroidal 2pi factor)
psiflu=permute(TRANSP.T.PLFLX,[2 1]); % on flux surface grid, [T m^2]
psiflu=psiflu(:,itransp);
nflxsurf=numel(psiflu); % number of flux surf.
[~,iflxsurf]=min(abs(psiflu-psinrm_in^2*psiflu(end))); % flux-grid point closest to user specified radial location

% Normalized sqrt(psi)
rho_TRANSP=zeros(1,nflxsurf);
for i=1:nflxsurf
    rho_TRANSP(i)=sqrt(TRANSP.T.PLFLX(itransp,i)/TRANSP.T.PLFLX(itransp,end));
end

% Radial coordinates in meters
a=(TRANSP.G.RMAJM(itransp,end)-TRANSP.G.RMAJM(itransp,1))/2.; % GS2 Lref
rho_GS2=zeros(1,nflxsurf); % GS2 definition of rho for irho=2, not yet normalized
for indx=1:nflxsurf
    rho_GS2(indx)=(TRANSP.G.RMAJM(itransp,nflxsurf+1+indx)-TRANSP.G.RMAJM(itransp,nflxsurf+1-indx))/2.;
end

%% Extract profiles from EFIT, EL and ION

% Normalized sqrt(psi) grid used in chain2
rho_chain2=linspace(0.,1.,21);

% Safety factor
q=interpol(EFIT.RMJO(ijp,:),EFIT.Q(ijp,:),TRANSP.G.RMAJM(itransp,nflxsurf+2:end));

% Temperature [eV]
ti=interpol(rho_chain2,ION.TI(ijp,:),rho_TRANSP);

% Tor. angular vel.
omega=interpol(rho_chain2,ION.ANGF(ijp,:),rho_TRANSP);
domega_drho=interpol(rho_GS2,omega,rho_GS2,1);

% Reference quantities
tref=ti(iflxsurf);
mref=2.*1.67e-27;
vtref=sqrt(2.*e*tref/mref);

% dist_fn_knobs
mach=omega(iflxsurf)*a/vtref;
g_exb=a/vtref*rho_GS2(iflxsurf)/q(iflxsurf)*domega_drho(iflxsurf);

end