function JETPEAK_to_gs2_input_old(ijp,psinrm_in,outfile_template,outfile_name,plot_verbose)

% loading databases
load databases/JETPEAK_2017_04.mat
load databases/TRANSP_2017_3.mat

if nargin < 5
    plot_verbose=0;
end

% Elementary charge
e=1.602e-19;

% Ignore NaN warnings
warning('off','MATLAB:chckxy:IgnoreNaN');

% TRANSP index corresponding to ijp
itransp=find(TRANSP.JPI==ijp);


%% Extract mag. geometry from TRANSP

% Poloidal flux (GS2 definition, i.e. without toroidal 2pi factor)
psi=permute(TRANSP.T.PSI,[2 3 1]); % on rectangular grid, [T m^2]
psi=psi(:,:,itransp);
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
Rmag=1.e-2*TRANSP.T.RAXIS(itransp);
Rflu=1.e-2*permute(TRANSP.T.RFLU,[2 3 1]); % flux surface grid
Rflu=Rflu(:,:,itransp);
Rpsi= 1.e-2*TRANSP.T.PSIR; % rectangular grid
a=(TRANSP.G.RMAJM(itransp,end)-TRANSP.G.RMAJM(itransp,1))/2.; % GS2 Lref
rho_GS2=zeros(1,nflxsurf); % GS2 definition of rho for irho=2, not yet normalized
for indx=1:nflxsurf
    rho_GS2(indx)=(TRANSP.G.RMAJM(itransp,nflxsurf+1+indx)-TRANSP.G.RMAJM(itransp,nflxsurf+1-indx))/2.;
end
Rmaj=TRANSP.G.RMAJM(itransp,nflxsurf+1+iflxsurf)-rho_GS2(iflxsurf);
Rgeo=Rmaj; % Free choice for Rgeo. This then defines Bref.
RMAJM=permute(TRANSP.G.RMAJM,[2 1]);
RMAJM=RMAJM(nflxsurf+2:end,itransp);
dR_drho=interpol(rho_GS2,RMAJM,rho_GS2,1);

% Vertical coord. in meters
Zmag=1.e-2*TRANSP.T.YAXIS(itransp);
Zflu=1.e-2*permute(TRANSP.T.ZFLU,[2 3 1]); 
Zflu=Zflu(:,:,itransp);
Zpsi= 1.e-2*TRANSP.T.PSIZ;

% Elongation
kappa=permute(TRANSP.T.ELONG,[2 1]);
kappa=kappa(:,itransp);
dkappa_drho=interpol(rho_GS2,kappa,rho_GS2,1);

% Triangularity
delta=permute(TRANSP.T.TRIANG,[2 1]);
delta=delta(:,itransp);
ddelta_drho=interpol(rho_GS2,delta,rho_GS2,1);


%% Extract profiles from EFIT, EL and ION

% Normalized sqrt(psi) grid used in chain2
rho_chain2=linspace(0.,1.,21);

% Safety factor
q=interpol(EFIT.RMJO(ijp,:),EFIT.Q(ijp,:),TRANSP.G.RMAJM(itransp,nflxsurf+2:end));
dq_drho=interpol(rho_GS2,q,rho_GS2,1);

% Density
ni=EL.NE(ijp,:)-6.0*ION.NC(ijp,:);
ni=interpol(rho_chain2,ni,rho_TRANSP);
dni_drho=interpol(rho_GS2,ni,rho_GS2,1);
ne=interpol(rho_chain2,EL.NE(ijp,:),rho_TRANSP);
dne_drho=interpol(rho_GS2,ne,rho_GS2,1);

% Temperature [eV]
ti=interpol(rho_chain2,ION.TI(ijp,:),rho_TRANSP);
dti_drho=interpol(rho_GS2,ti,rho_GS2,1);
te=interpol(rho_chain2,EL.TE(ijp,:),rho_TRANSP);
dte_drho=interpol(rho_GS2,te,rho_GS2,1);

% Plasma pressure
p=e*ni.*ti+e*ne.*te;
dp_drho=interpol(rho_GS2,p,rho_GS2,1);

% Tor. angular vel.
omega=interpol(rho_chain2,ION.ANGF(ijp,:),rho_TRANSP);
domega_drho=interpol(rho_GS2,omega,rho_GS2,1);

% Effective charge
warning('PLEASE CHECK DEFINITION OF ZEFF !')
Zeff=BASIC.ZEFV(ijp);


%% Calculate I(psi) to determine Bref for given Rgeo
% where I(psi) comes from expressing B as B = I*grad(zeta) + grad(zeta) x grad(psi)

Rflu_loc=permute(Rflu,[2 1]); 
Rflu_loc=Rflu_loc(:,iflxsurf); % Rflu for all poloidal locations on chosen flux surf
Zflu_loc=permute(Zflu,[2 1]);
Zflu_loc=Zflu_loc(:,iflxsurf); % Zflu for all poloidal locations on chosen flux surf
I=compute_I(Rmag,Zmag,Rflu_loc,Zflu_loc,Rpsi,Zpsi,psi,q(iflxsurf));
% In GS2, Bref is defined via Bref=I/(Rgeo*a)
Bref=I/Rgeo; % no 1/a since Rgeo not normalised yet


%% Collect into jet_out structure

jet_out.a=a;
jet_out.Rmaj=Rmaj;
jet_out.Rgeo=Rgeo;
jet_out.rho=rho_GS2(iflxsurf);
jet_out.dR_drho=dR_drho(iflxsurf);
jet_out.q=q(iflxsurf);
jet_out.dq_drho=dq_drho(iflxsurf);
jet_out.kappa=kappa(iflxsurf);
jet_out.dkappa_drho=dkappa_drho(iflxsurf);
jet_out.delta=delta(iflxsurf);
jet_out.ddelta_drho=ddelta_drho(iflxsurf);
jet_out.Bref=Bref;
jet_out.ni=ni(iflxsurf);
jet_out.dni_drho=dni_drho(iflxsurf);
jet_out.ne=ne(iflxsurf);
jet_out.dne_drho=dne_drho(iflxsurf);
jet_out.ti=ti(iflxsurf);
jet_out.dti_drho=dti_drho(iflxsurf);
jet_out.te=te(iflxsurf);
jet_out.dte_drho=dte_drho(iflxsurf);
jet_out.p=p(iflxsurf);
jet_out.dp_drho=dp_drho(iflxsurf);
jet_out.omega=omega(iflxsurf);
jet_out.domega_drho=domega_drho(iflxsurf);
jet_out.Zeff=Zeff;


%% GS2-normalisation

gs2_in=normalise_jet_out(jet_out);


%% Generate input file

generate_gs2_infile(gs2_in,outfile_template,outfile_name);


%% Check data prof., interpolations, derivatives (if plot_verbose)

if plot_verbose
    
    % Major radius
    figure;
    plot(rho_TRANSP,RMAJM,'bx')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('R')
    figure;
    plot(rho_TRANSP,dR_drho,'m-x')
    legend('GS2 d/drho')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$\frac{dR}{d\rho}$','Interpreter','LaTex')
    
    % Safety factor
    figure;
    plot(EFIT.RMJO(ijp,:),EFIT.Q(ijp,:),'bx')
    hold on
    plot(TRANSP.G.RMAJM(itransp,nflxsurf+2:end),q,'r-x')
    legend('JETPEAK','GS2 fit')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$q$','Interpreter','LaTex')
    figure;
    plot(rho_TRANSP,dq_drho,'m-x')
    legend('GS2 d/drho')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$\frac{dq}{d\rho}$','Interpreter','LaTex')
    
    % Elongation
    figure;
    plot(rho_TRANSP,TRANSP.T.ELONG(itransp,:),'bx')
    hold on
    plot(rho_TRANSP,kappa,'r-x')    
    legend('JETPEAK','GS2 fit')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$\kappa$','Interpreter','LaTex')
    figure;
    plot(rho_TRANSP,dkappa_drho,'m-x')
    legend('GS2 d/drho')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$\frac{d\kappa}{d\rho}$','Interpreter','LaTex')
    
    % Triangularity
    figure;
    plot(rho_TRANSP,TRANSP.T.TRIANG(itransp,:),'bx')
    hold on
    plot(rho_TRANSP,delta,'r-x')    
    legend('JETPEAK','GS2 fit')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$\delta$','Interpreter','LaTex')
    figure;
    plot(rho_TRANSP,ddelta_drho,'m-x')
    legend('GS2 d/drho')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$\frac{d\delta}{d\rho}$','Interpreter','LaTex')
    
    % Ion density
    figure;
    plot(rho_chain2,EL.NE(ijp,:)-6.0*ION.NC(ijp,:),'bx')
    hold on
    plot(rho_TRANSP,ni,'r-x')    
    legend('JETPEAK','GS2 fit')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$n_i$','Interpreter','LaTex')
    figure;
    plot(rho_TRANSP,dni_drho,'m-x')
    legend('GS2 d/drho')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$\frac{dn_i}{d\rho}$','Interpreter','LaTex')
    
    % Electron density
    figure;
    plot(rho_chain2,EL.NE(ijp,:),'bx')
    hold on
    plot(rho_TRANSP,ne,'r-x')    
    legend('JETPEAK','GS2 fit')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$n_e$','Interpreter','LaTex')
    figure;
    plot(rho_TRANSP,dne_drho,'m-x')
    legend('GS2 d/drho')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$\frac{dn_e}{d\rho}$','Interpreter','LaTex')
    
    % Ion temperature
    figure;
    plot(rho_chain2,ION.TI(ijp,:),'bx')
    hold on
    plot(rho_TRANSP,ti,'r-x')    
    legend('JETPEAK','GS2 fit')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$T_i$','Interpreter','LaTex')
    figure;
    plot(rho_TRANSP,dti_drho,'m-x')
    legend('GS2 d/drho')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$\frac{dT_i}{d\rho}$','Interpreter','LaTex')
    
    % Electron temperature
    figure;
    plot(rho_chain2,EL.TE(ijp,:),'bx')
    hold on
    plot(rho_TRANSP,te,'r-x')    
    legend('JETPEAK','GS2 fit')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$T_e$','Interpreter','LaTex')
    figure;
    plot(rho_TRANSP,dte_drho,'m-x')
    legend('GS2 d/drho')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$\frac{dT_e}{d\rho}$','Interpreter','LaTex')
    
    % Toroidal angular velocity
    figure;
    plot(rho_chain2,ION.ANGF(ijp,:),'bx')
    hold on
    plot(rho_TRANSP,omega,'r-x')    
    legend('JETPEAK','GS2 fit')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$\omega$','Interpreter','LaTex')
    figure;
    plot(rho_TRANSP,domega_drho,'m-x')
    legend('GS2 d/drho')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$\frac{d\omega}{d\rho}$','Interpreter','LaTex')
    
end


end

