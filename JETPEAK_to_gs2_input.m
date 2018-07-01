function JETPEAK_to_gs2_input(ijp,psinrm_in,outfile_template,outfile_name,plot_verbose)

% loading databases
load ~/Dropbox/oxford/physics/jetpeak/databases/JETPEAK_2017_04_1661torq.mat
load ~/Dropbox/oxford/physics/jetpeak/databases/TRANSP_2017_3.mat

warning('ONLY HAVE NB TORQUE DEPOSITION FOR ijp=1661 ...')

if nargin < 5
    plot_verbose=0;
end

% Elementary charge
e=1.602e-19;

% Ignore NaN warnings
warning('off','MATLAB:chckxy:IgnoreNaN');

% TRANSP index corresponding to ijp
itransp=find(TRANSP.JPI==ijp);


%% Choosing gs2 zeta to be co-current
% see notes about signs
sign_bpol = 1; % choice we make, means zeta_gs2 is in Ip direction
sign_btor = sign(BASIC.BT(ijp)*BASIC.IP(ijp));
sign_q = sign_bpol*sign_btor;
sign_Rgeo = sign_q;


%% Extract mag. geometry from TRANSP

% Poloidal flux (GS2 definition, i.e. without toroidal 2pi factor)
psi=permute(TRANSP.T.PSI,[2 3 1]); % on rectangular grid, [T m^2]
psi=psi(:,:,itransp);
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
Rmag=1.e-2*TRANSP.T.RAXIS(itransp);
Rflu=1.e-2*permute(TRANSP.T.RFLU,[2 3 1]); % flux surface grid
Rflu=Rflu(:,:,itransp);
Rpsi= 1.e-2*TRANSP.T.PSIR; % rectangular grid
a=(TRANSP.G.RMAJM(itransp,end)-TRANSP.G.RMAJM(itransp,1))/2.; % GS2 Lref
rhoc_TRANSP=zeros(1,nflxsurf); % GS2 definition of rho for irho=2, not yet normalized
Rmaj=zeros(1,nflxsurf); % Rmaj definition for iflux ~= 1, not yet normalized
for indx=1:nflxsurf
    rhoc_TRANSP(indx)=(TRANSP.G.RMAJM(itransp,nflxsurf+1+indx)-TRANSP.G.RMAJM(itransp,nflxsurf+1-indx))/2.;
    Rmaj(indx)=TRANSP.G.RMAJM(itransp,nflxsurf+1+indx)-rhoc_TRANSP(indx);
end
Rgeo=sign_Rgeo*Rmaj(iflxsurf); % Free choice for Rgeo. This then defines Bref., see notes about signs
dR_drho=interpol(rhoc_TRANSP,Rmaj,rhoc_TRANSP,1);

% Vertical coord. in meters
Zmag=1.e-2*TRANSP.T.YAXIS(itransp);
Zflu=1.e-2*permute(TRANSP.T.ZFLU,[2 3 1]);
Zflu=Zflu(:,:,itransp);
Zpsi= 1.e-2*TRANSP.T.PSIZ;

% Elongation
kappa=permute(TRANSP.T.ELONG,[2 1]);
kappa=kappa(:,itransp);
dkappa_drho=interpol(rhoc_TRANSP,kappa,rhoc_TRANSP,1);

% Triangularity
delta=permute(TRANSP.T.TRIANG,[2 1]);
delta=asin(delta(:,itransp));
ddelta_drho=interpol(rhoc_TRANSP,delta,rhoc_TRANSP,1);


%% Extract profiles from EFIT, EL and ION

% Normalized sqrt(psi) grid used in chain2
sqrt_psin_chain2=linspace(0.,1.,21);

% Safety factor with adjusted
q = sign_q * interpol(EFIT.RMJO(ijp,:),abs(EFIT.Q(ijp,:)),TRANSP.G.RMAJM(itransp,nflxsurf+2:end));
dq_drho=interpol(rhoc_TRANSP,q,rhoc_TRANSP,1);

% Density
ni=EL.NE(ijp,:)-6.0*ION.NC(ijp,:);
ni=interpol(sqrt_psin_chain2,ni,sqrt_psin_TRANSP);
dni_drho=interpol(rhoc_TRANSP,ni,rhoc_TRANSP,1);
ne=interpol(sqrt_psin_chain2,EL.NE(ijp,:),sqrt_psin_TRANSP);
dne_drho=interpol(rhoc_TRANSP,ne,rhoc_TRANSP,1);

% Temperature [eV]
ti=interpol(sqrt_psin_chain2,ION.TI(ijp,:),sqrt_psin_TRANSP);
dti_drho=interpol(rhoc_TRANSP,ti,rhoc_TRANSP,1);
te=interpol(sqrt_psin_chain2,EL.TE(ijp,:),sqrt_psin_TRANSP);
dte_drho=interpol(rhoc_TRANSP,te,rhoc_TRANSP,1);

% Plasma pressure
p=e*ni.*ti+e*ne.*te;
dp_drho=interpol(rhoc_TRANSP,p,rhoc_TRANSP,1);

% Tor. angular vel.
omega = interpol(sqrt_psin_chain2,ION.ANGF(ijp,:),sqrt_psin_TRANSP);
domega_drho = interpol(rhoc_TRANSP,omega,rhoc_TRANSP,1);
sign_omega = sign(omega(iflxsurf))*sign(BASIC.IP(ijp)); % see notes about signs
sign_domega_drho = sign_q*sign(domega_drho(iflxsurf))*sign(BASIC.IP(ijp)); % see notes about signs
omega = sign_omega*omega;
domega_drho = sign_domega_drho*domega_drho;

% Effective charge
% Assume only impurity species is Carbon
Zimp=6.;
Zeff=1.+(Zimp^2-Zimp)*ION.NC(ijp,:)./EL.NE(ijp,:);
Zeff=interpol(sqrt_psin_chain2,Zeff,sqrt_psin_TRANSP);


%% Calculate I(psi) to determine Bref for given Rgeo
% where I(psi) comes from expressing B as B = I*grad(zeta) + grad(zeta) x grad(psi)

Rflu_loc=permute(Rflu,[2 1]); 
Rflu_loc=Rflu_loc(:,iflxsurf); % Rflu for all poloidal locations on chosen flux surf
Zflu_loc=permute(Zflu,[2 1]);
Zflu_loc=Zflu_loc(:,iflxsurf); % Zflu for all poloidal locations on chosen flux surf
I = compute_I(Rmag,Zmag,Rflu_loc,Zflu_loc,Rpsi,Zpsi,psi,q(iflxsurf));
% In GS2, Bref is defined via Bref=I/(Rgeo*a)
Bref=abs(I/Rgeo); % no 1/a since Rgeo not normalised yet, see notes about signs


%% Collect into jet_out structure

jet_out.a=a;
jet_out.nref=ni(iflxsurf); % reference dens. and temp. from ions
jet_out.tref=ti(iflxsurf);
jet_out.mref=1.675e-27+1.6726e-27;
jet_out.Bref=Bref;

jet_out.Rmaj=Rmaj(iflxsurf);
jet_out.Rgeo=Rgeo;
jet_out.rho=rhoc_TRANSP(iflxsurf);
jet_out.dR_drho=dR_drho(iflxsurf);
jet_out.q=q(iflxsurf);
jet_out.dq_drho=dq_drho(iflxsurf);
jet_out.kappa=kappa(iflxsurf);
jet_out.dkappa_drho=dkappa_drho(iflxsurf);
jet_out.delta=delta(iflxsurf);
jet_out.ddelta_drho=ddelta_drho(iflxsurf);
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
jet_out.Zeff=Zeff(iflxsurf);


%% GS2-normalisation

gs2_in=normalise_jet_out(jet_out);


%% Generate input file

generate_gs2_infile(gs2_in,outfile_template,outfile_name);


%% Check data prof., interpolations, derivatives (if plot_verbose)

if plot_verbose
    
    % Major radius
    figure;
    plot(sqrt_psin_TRANSP,Rmaj,'bx')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('R')
    figure;
    plot(sqrt_psin_TRANSP,dR_drho,'m-x')
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
    plot(sqrt_psin_TRANSP,dq_drho,'m-x')
    legend('GS2 d/drho')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$\frac{dq}{d\rho}$','Interpreter','LaTex')
    
    % Elongation
    figure;
    plot(sqrt_psin_TRANSP,TRANSP.T.ELONG(itransp,:),'bx')
    hold on
    plot(sqrt_psin_TRANSP,kappa,'r-x')    
    legend('JETPEAK','GS2 fit')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$\kappa$','Interpreter','LaTex')
    figure;
    plot(sqrt_psin_TRANSP,dkappa_drho,'m-x')
    legend('GS2 d/drho')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$\frac{d\kappa}{d\rho}$','Interpreter','LaTex')
    
    % Triangularity
    figure;
    plot(sqrt_psin_TRANSP,TRANSP.T.TRIANG(itransp,:),'bx')
    hold on
    plot(sqrt_psin_TRANSP,delta,'r-x')    
    legend('JETPEAK','GS2 fit')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$\delta$','Interpreter','LaTex')
    figure;
    plot(sqrt_psin_TRANSP,ddelta_drho,'m-x')
    legend('GS2 d/drho')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$\frac{d\delta}{d\rho}$','Interpreter','LaTex')
    
    % Ion density
    figure;
    plot(sqrt_psin_chain2,EL.NE(ijp,:)-6.0*ION.NC(ijp,:),'bx')
    hold on
    plot(sqrt_psin_TRANSP,ni,'r-x')    
    legend('JETPEAK','GS2 fit')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$n_i$','Interpreter','LaTex')
    figure;
    plot(sqrt_psin_TRANSP,dni_drho,'m-x')
    legend('GS2 d/drho')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$\frac{dn_i}{d\rho}$','Interpreter','LaTex')
    
    % Electron density
    figure;
    plot(sqrt_psin_chain2,EL.NE(ijp,:),'bx')
    hold on
    plot(sqrt_psin_TRANSP,ne,'r-x')    
    legend('JETPEAK','GS2 fit')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$n_e$','Interpreter','LaTex')
    figure;
    plot(sqrt_psin_TRANSP,dne_drho,'m-x')
    legend('GS2 d/drho')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$\frac{dn_e}{d\rho}$','Interpreter','LaTex')
    
    % Ion temperature
    figure;
    plot(sqrt_psin_chain2,ION.TI(ijp,:),'bx')
    hold on
    plot(sqrt_psin_TRANSP,ti,'r-x')    
    legend('JETPEAK','GS2 fit')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$T_i$','Interpreter','LaTex')
    figure;
    plot(sqrt_psin_TRANSP,dti_drho,'m-x')
    legend('GS2 d/drho')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$\frac{dT_i}{d\rho}$','Interpreter','LaTex')
    
    % Electron temperature
    figure;
    plot(sqrt_psin_chain2,EL.TE(ijp,:),'bx')
    hold on
    plot(sqrt_psin_TRANSP,te,'r-x')    
    legend('JETPEAK','GS2 fit')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$T_e$','Interpreter','LaTex')
    figure;
    plot(sqrt_psin_TRANSP,dte_drho,'m-x')
    legend('GS2 d/drho')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$\frac{dT_e}{d\rho}$','Interpreter','LaTex')
    
    % Toroidal angular velocity
    figure;
    plot(sqrt_psin_chain2,ION.ANGF(ijp,:),'bx')
    hold on
    plot(sqrt_psin_TRANSP,omega,'r-x')    
    legend('JETPEAK','GS2 fit')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$\omega$','Interpreter','LaTex')
    figure;
    plot(sqrt_psin_TRANSP,domega_drho,'m-x')
    legend('GS2 d/drho')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$\frac{d\omega}{d\rho}$','Interpreter','LaTex')
    
end


end

