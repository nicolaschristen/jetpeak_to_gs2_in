load databases/JETPEAK_2017_04.mat
load databases/TRANSP_2017_3.mat

% Ignore NaN warnings
warning('off','MATLAB:chckxy:IgnoreNaN');

ijp_plots=[];
mach_edge=[];
mach_in=[];
g_exb_edge=[];
g_exb_in=[];
for ijp=1:1847
    itransp=find(TRANSP.JPI==ijp);
    if ~isempty(itransp) && isempty(find(isnan(EFIT.RMJO(ijp,:)),1)) ...
            && isempty(find(isnan(ION.ANGF(ijp,:)),1))
        itransp=itransp(1);
        
        if isempty(find(isnan(TRANSP.G.RMAJM(itransp,:)),1)) ...
                && isempty(find(isnan(TRANSP.T.PLFLX(itransp,:)),1))
            
            ijp_plots=[ijp_plots,ijp];
            
            % compute GS2 mach & g_exb at the edge border
            psinrm_edge=0.8;
            [m_edge,g_edge]=compute_mach_g_exb(ijp,itransp,psinrm_edge,TRANSP,EFIT,ION);
            mach_edge=[mach_edge,m_edge];
            g_exb_edge=[g_exb_edge,g_edge];
            
            % same for inner border
            psinrm_in=0.3;
            [m_in,g_in]=compute_mach_g_exb(ijp,itransp,psinrm_in,TRANSP,EFIT,ION);
            mach_in=[mach_in,m_in];
            g_exb_in=[g_exb_in,g_in];
            
        end
              
    end    
end

figure
plot(ijp_plots,mach_edge,'b*')
xlabel('ijp')
ylabel('mach [vthref/a]')
title('GS2 $\omega$ at $\sqrt{\psi_N}=0.8$','Interpreter','Latex')

figure
plot(ijp_plots,g_exb_edge,'b*')
xlabel('ijp')
ylabel('g\_exb [q*vthref/(a*rho)]')
title('GS2 $\frac{d\omega}{d\rho}$ at $\sqrt{\psi_N}=0.8$','Interpreter','Latex')

figure
plot(ijp_plots,mach_in,'b*')
xlabel('ijp')
ylabel('mach [vthref/a]')
title('GS2 $\omega$ at $\sqrt{\psi_N}=0.3$','Interpreter','Latex')

figure
plot(ijp_plots,g_exb_in,'b*')
xlabel('ijp')
ylabel('g\_exb [q*vthref/(a*rho)]')
title('GS2 $\frac{d\omega}{d\rho}$ at $\sqrt{\psi_N}=0.3$','Interpreter','Latex')

function [mach,g_exb] = compute_mach_g_exb(ijp,itransp,psinrm_in,TRANSP,EFIT,ION)

% Elementary charge
e=1.602e-19;

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