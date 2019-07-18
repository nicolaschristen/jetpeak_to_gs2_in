%% Plot mach and g_exb to select interesting shot from JETPEAK
%
% This script plots mach & g_exb at sqrt(psipol/psipol_LCFS)=0.3 and 0.8,
% with:
%       mach = Omega_phi * rref/vtref
%       g_exb = dOmega_phi/dr_psi * r_psi/q * rref/vtref

load ../databases/JETPEAK_2017_04_1661torq.mat
load ../databases/TRANSP_2017_3.mat

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