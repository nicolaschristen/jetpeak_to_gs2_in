%% Generate GS2 input files from combined JETPEAK and TRANSP databases
%
% Input :   ijp -- shot index in JETPEAK DB
%           psinrm_in -- sqrt(psipol/psipol_LCFS) of selected flux surf.
%           infile_template -- GS2 input file template
%           new_infile_name -- name of input file to be generated
%           plot_verbose_main (=0) -- create plots to check params & fits
%           plot_verbose_Ipsi_integration (=0)
%                      -- create plots to check fits in I(psi) integration
%
% Output:   -
%
% Plasma parameters are evaluated at the flux surface
% nearest to psinrm_in (no interpolation).
%
% All quantities are in SI units,
% except for temperatures which are in eV.
%
% Plots can be generated to visually check parameters and fits:
% below, set plot_verbose_main = 1,
% and/or plot_verbose_Ipsi_integration = 0
%
function JETPEAK_to_gs2_input(ijp,psinrm_in,infile_template,new_infile_name, ...
    plot_verbose_main, plot_verbose_Ipsi_integration)

if nargin < 6
    % Set this to 1 to check fits in computation of I(psi)
    plot_verbose_Ipsi_integration = 0;
    if nargin < 5
        % Set this to 1 to check fits of plasma parameters
        plot_verbose_main = 0;
    end
end

% loading databases
load databases/JETPEAK_2017_04_1661torq.mat
load databases/TRANSP_2017_3.mat

fprintf('\nOnly have NB torque deposition for ijp=1661. \n\n')

% Ignore NaN warnings
warning('off','MATLAB:chckxy:IgnoreNaN');

% TRANSP index corresponding to ijp
itransp=find(TRANSP.JPI==ijp);

% Physical constants
e = 1.602e-19; % elementary charge
eps0 = 8.8541878128e-12; % vacuum permittivity
mp = 1.673e-27; % proton mass
me = 9.109e-31; % electron mass

% Assume Deuterium plasma with Carbon impurity
Zmain=1.;
Zimp=6.;


%% Signs (see notes)
sign_bpol = 1; % choice we make, means zeta_gs2 is in Ip direction
sign_btor = sign(BASIC.BT(ijp)*BASIC.IP(ijp));
sign_q = sign_bpol*sign_btor;
sign_Rgeo = sign_q;


%% Extract mag. geometry from TRANSP

% Poloidal flux [T m^2]
% with definition equivalent to GS2's, i.e. without toroidal 2pi factor.
% See Edmund's thesis and the TRANSP variable list.
psi=permute(TRANSP.T.PSI,[3 2 1]); % on rectangular grid
psi=psi(:,:,itransp); % psi(iR,iZ)
psiflu=permute(TRANSP.T.PLFLX,[2 1]); % on flux surface grid
psiflu=psiflu(:,itransp); % psiflu(iflxsurf)
nflxsurf=numel(psiflu); % number of flux surf.
% flux-grid point closest to user specified radial location
[~,iflxsurf]=min(abs(psiflu-psinrm_in^2*psiflu(end)));

% Normalized sqrt(psi) []
sqrt_psin_TRANSP=zeros(1,nflxsurf);
for i=1:nflxsurf
    sqrt_psin_TRANSP(i)=sqrt(TRANSP.T.PLFLX(itransp,i)/TRANSP.T.PLFLX(itransp,end));
end

% Radial coordinates [m]
Rmag=1.e-2*TRANSP.T.RAXIS(itransp);
Rflu=1.e-2*permute(TRANSP.T.RFLU,[2 3 1]); % flux surface grid
Rflu=Rflu(:,:,itransp); % (iflx,itheta)
Rpsi= 1.e-2*TRANSP.T.PSIR; % rectangular grid
a=(TRANSP.G.RMAJM(itransp,end)-TRANSP.G.RMAJM(itransp,1))/2.; % GS2 Lref
rpsi_TRANSP=zeros(1,nflxsurf); % GS2 definition of rho for irho=2, not yet normalized
Rmaj=zeros(1,nflxsurf); % Rmaj definition for iflux ~= 1, not yet normalized
for indx=1:nflxsurf
    rpsi_TRANSP(indx)= ...
        (TRANSP.G.RMAJM(itransp,nflxsurf+1+indx) ...
           - TRANSP.G.RMAJM(itransp,nflxsurf+1-indx))/2.;
    Rmaj(indx)=TRANSP.G.RMAJM(itransp,nflxsurf+1+indx)-rpsi_TRANSP(indx);
end
% Free choice for Rgeo. This then defines Bref., see notes about signs
Rgeo=sign_Rgeo*abs(Rmaj(iflxsurf));
dR_drho=interpol(rpsi_TRANSP,Rmaj,rpsi_TRANSP,1);

% Vertical coord. [m]
Zmag=1.e-2*TRANSP.T.YAXIS(itransp);
Zflu=1.e-2*permute(TRANSP.T.ZFLU,[2 3 1]);
Zflu=Zflu(:,:,itransp); % (iflx,itheta)
Zpsi= 1.e-2*TRANSP.T.PSIZ;

% Elongation []
kappa=permute(TRANSP.T.ELONG,[2 1]);
kappa=kappa(:,itransp);
dkappa_drho=interpol(rpsi_TRANSP,kappa,rpsi_TRANSP,1);

% Triangularity []
delta=permute(TRANSP.T.TRIANG,[2 1]);
delta=asin(delta(:,itransp));
ddelta_drho=interpol(rpsi_TRANSP,delta,rpsi_TRANSP,1);


%% Extract q and plasma profiles from EFIT, EL and ION

% Normalized sqrt(psi) grid used in chain2 []
sqrt_psin_chain2=linspace(0.,1.,21);

% Safety factor []
q = sign_q * interpol(EFIT.RMJO(ijp,:), ...
    abs(EFIT.Q(ijp,:)),TRANSP.G.RMAJM(itransp,nflxsurf+2:end));
dq_drho=interpol(rpsi_TRANSP,q,rpsi_TRANSP,1);

% Density [m^{-3}]
ni=EL.NE(ijp,:)-6.0*ION.NC(ijp,:);
ni=interpol(sqrt_psin_chain2,ni,sqrt_psin_TRANSP);
nc=interpol(sqrt_psin_chain2,ION.NC(ijp,:),sqrt_psin_TRANSP);
ne=interpol(sqrt_psin_chain2,EL.NE(ijp,:),sqrt_psin_TRANSP);
fprintf('Should Carbon be included ? nc/ni = %1.4f\n',nc(iflxsurf)/ni(iflxsurf))
if nc(iflxsurf)/ni(iflxsurf) >= 0.05
    fprintf('\tnc/ni >= 0.05: including carbon.\n\n')
    add_carbon = 1;
else
    fprintf('\tnc/ni < 0.05: not including carbon.\n')
    fprintf('\tni is adapted so that QN holds.\n\n')
    add_carbon = 0;
    ni = ni + Zimp/Zmain*nc;
end
dni_drho=interpol(rpsi_TRANSP,ni,rpsi_TRANSP,1);
dne_drho=interpol(rpsi_TRANSP,ne,rpsi_TRANSP,1);
dnc_drho=interpol(rpsi_TRANSP,nc,rpsi_TRANSP,1);

% Temperature [eV]
ti=interpol(sqrt_psin_chain2,ION.TI(ijp,:),sqrt_psin_TRANSP);
dti_drho=interpol(rpsi_TRANSP,ti,rpsi_TRANSP,1);
te=interpol(sqrt_psin_chain2,EL.TE(ijp,:),sqrt_psin_TRANSP);
dte_drho=interpol(rpsi_TRANSP,te,rpsi_TRANSP,1);
if add_carbon
    fprintf('T_C is set to T_i since there are no direct measurements.\n\n')
    tc=ti; % no direct measurement for T_imp
    dtc_drho=dti_drho;
end

% Plasma pressure [Pa]
p=e*ni.*ti+e*ne.*te;
if add_carbon
    p = p+e*nc.*tc;
end
dp_drho=interpol(rpsi_TRANSP,p,rpsi_TRANSP,1);

% Collisionality (See notes) [s^{-1}]
colfac = e^(-3/2)/(4*pi*eps0)^2; % conversion factor from cgs to SI+eV
Zi = 1; % main ion is deuterium
Ze = -1;
mi = 2*mp; % main ion is deuterium
loglamda = coulomb_log(Zi,mi,ne,ni,te,ti);
nu_ii = colfac * 4*pi*Zi^4*e^4*ni.*loglamda./(sqrt(mi)*(2*ti).^(3/2));
nu_ee = colfac * 4*pi*Ze^4*e^4*ne.*loglamda./(sqrt(me)*(2*te).^(3/2));
if add_carbon
    Zc = 6;
    mc = 12*mp;
    loglamda_c = coulomb_log(Zc,mc,ne,nc,te,tc);
    nu_cc = colfac * 4*pi*Zc^4*e^4*nc.*loglamda_c./(sqrt(mc)*(2*tc).^(3/2));
end

% Tor. angular vel. [radian/s]
omega = interpol(sqrt_psin_chain2,ION.ANGF(ijp,:),sqrt_psin_TRANSP);
domega_drho = interpol(rpsi_TRANSP,omega,rpsi_TRANSP,1);
% see notes about signs
sign_mach = sign(omega(iflxsurf))*sign(BASIC.IP(ijp));
% see notes about signs
sign_domega_drho = sign(domega_drho(iflxsurf))*sign(BASIC.IP(ijp));
omega = sign_mach*abs(omega);
domega_drho = sign_domega_drho*abs(domega_drho);

% Effective charge [C]
if add_carbon
   Zeff = (Zmain^2*ni+Zimp^2*nc)./(Zmain*ni+Zimp*nc);
else
   Zeff = Zmain*ones(nflxsurf,1);
end


%% Calculate I(psi) to determine Bref for given Rgeo
% where I(psi) comes from expressing B as B = I*grad(zeta) + grad(zeta) x grad(psi)

Rflu_loc=permute(Rflu,[2 1]); 
Rflu_loc=Rflu_loc(:,iflxsurf); % Rflu for all poloidal locations on chosen flux surf
Zflu_loc=permute(Zflu,[2 1]);
Zflu_loc=Zflu_loc(:,iflxsurf); % Zflu for all poloidal locations on chosen flux surf
I = compute_I(Rmag,Zmag,Rflu_loc,Zflu_loc,Rpsi,Zpsi,psi,q(iflxsurf), ...
    psiflu(iflxsurf), plot_verbose_Ipsi_integration);
% In GS2, Bref is defined via Bref=I/(Rgeo*a)
Bref=abs(I/Rgeo); % no 1/a since Rgeo not normalised yet, see notes about signs


%% Collect into jet_out structure

jet_out.a=a;
jet_out.nref=ni(iflxsurf); % reference dens. and temp. from ions
jet_out.tref=ti(iflxsurf);
jet_out.mref=1.675e-27+1.6726e-27; % assume main ion is deuterium
jet_out.Bref=Bref;

jet_out.Rmaj=Rmaj(iflxsurf);
jet_out.Rgeo=Rgeo;
jet_out.rho=rpsi_TRANSP(iflxsurf);
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
if add_carbon
    jet_out.nc=nc(iflxsurf);
    jet_out.dnc_drho=dnc_drho(iflxsurf);
end
jet_out.ti=ti(iflxsurf);
jet_out.dti_drho=dti_drho(iflxsurf);
jet_out.te=te(iflxsurf);
jet_out.dte_drho=dte_drho(iflxsurf);
if add_carbon
    jet_out.tc=tc(iflxsurf);
    jet_out.dtc_drho=dtc_drho(iflxsurf);
end
jet_out.p=p(iflxsurf);
jet_out.dp_drho=dp_drho(iflxsurf);
jet_out.nu_ii=nu_ii(iflxsurf);
jet_out.nu_ee=nu_ee(iflxsurf);
if add_carbon
    jet_out.nu_cc=nu_cc(iflxsurf);
end
jet_out.omega=omega(iflxsurf);
jet_out.domega_drho=domega_drho(iflxsurf);
jet_out.Zeff=Zeff(iflxsurf);


%% GS2-normalisation

gs2_in=normalise_jet_out(jet_out, add_carbon);


%% Generate input file

generate_gs2_infile(gs2_in,infile_template,new_infile_name, add_carbon);


%% Check data prof., interpolations, derivatives (if plot_verbose)

if plot_verbose_main
    
    % Major radius
    figure;
    plot(rpsi_TRANSP,Rmaj,'b-x')
    hold on
    plot(rpsi_TRANSP, ...
        Rmaj(iflxsurf)+ ...
           (rpsi_TRANSP-rpsi_TRANSP(iflxsurf))*dR_drho(iflxsurf),'r-')
    legend('Experimental data','Tangent at point of interest')
    xlabel('$r_\psi$','Interpreter','LaTex')
    ylabel('$R_\psi$','Interpreter','LaTex')
    hold on
    plot(rpsi_TRANSP(iflxsurf),Rmaj(iflxsurf),'ro')
    hold off
    
    % Elongation
    figure;
    plot(rpsi_TRANSP,kappa,'b-x')
    hold on
    plot(rpsi_TRANSP, ...
        kappa(iflxsurf)+ ...
           (rpsi_TRANSP-rpsi_TRANSP(iflxsurf))*dkappa_drho(iflxsurf),'r-')
    legend('Experimental data','Tangent at point of interest')
    xlabel('$r_\psi$','Interpreter','LaTex')
    ylabel('$\kappa$','Interpreter','LaTex')
    hold on
    plot(rpsi_TRANSP(iflxsurf),kappa(iflxsurf),'ro')
    hold off
    
    % Triangularity
    figure;
    plot(rpsi_TRANSP,delta,'b-x')
    hold on
    plot(rpsi_TRANSP, ...
        delta(iflxsurf)+...
           (rpsi_TRANSP-rpsi_TRANSP(iflxsurf))*ddelta_drho(iflxsurf),'r-')
    legend('Experimental data','Tangent at point of interest')
    xlabel('$r_\psi$','Interpreter','LaTex')
    ylabel('$\delta$','Interpreter','LaTex')
    hold on
    plot(rpsi_TRANSP(iflxsurf),delta(iflxsurf),'ro')
    hold off
    % Check if delta_GS2 = asin(delta_TRANSP)
    figure;
    plot(rpsi_TRANSP,TRANSP.T.TRIANG(itransp,:),'b-x')
    hold on
    % Calculate with Ztop
    [~,ithet_top] = max(Zflu,[],2);
    R_top = zeros(nflxsurf,1);
    for iflx = 1:nflxsurf
        R_top(iflx) = Rflu(iflx,ithet_top(iflx));
    end
    R_top = R_top.';
    plot(rpsi_TRANSP,(Rmaj-R_top)./rpsi_TRANSP,'r-o')
    hold on
    % Calculate with Zbot
    [~,ithet_bot] = min(Zflu,[],2);
    R_bot = zeros(nflxsurf,1);
    for iflx = 1:nflxsurf
        R_bot(iflx) = Rflu(iflx,ithet_bot(iflx));
    end
    R_bot = R_bot.';
    plot(rpsi_TRANSP,(Rmaj-R_bot)./rpsi_TRANSP,'g-o')
    lgd = legend('$\delta_{TRANSP}$',...
        '$(R_\psi-R_{top})/r_\psi$',...
        '$(R_\psi-R_{bot})/r_\psi$');
    lgd.FontSize = 18;
    lgd.Location = 'northwest';
    xlabel('$r_\psi$','Interpreter','LaTex')
    hold off
    
    % Safety factor
    % Check q interpolation to TRANSP grid
    figure;
    plot(EFIT.RMJO(ijp,:),abs(EFIT.Q(ijp,:)),'bx')
    hold on
    plot(TRANSP.G.RMAJM(itransp,nflxsurf+2:end),abs(q),'r-')
    legend('EFIT','fit to TRANSP')
    xlabel('$R$','Interpreter','LaTex')
    ylabel('$\vert q\vert$','Interpreter','LaTex')
    title(['sgn(q)= ' num2str(sign_q)], 'FontSize',18)
    % Check derivative
    figure;
    plot(rpsi_TRANSP,q,'b-x')
    hold on
    plot(rpsi_TRANSP,...
        q(iflxsurf)+...
           (rpsi_TRANSP-rpsi_TRANSP(iflxsurf))*dq_drho(iflxsurf),'r-')
    legend('Experimental data','Tangent at point of interest')
    xlabel('$r_\psi$','Interpreter','LaTex')
    ylabel('$q$','Interpreter','LaTex')
    hold on
    plot(rpsi_TRANSP(iflxsurf),q(iflxsurf),'ro')
    hold off
    
    % Ion density
    % First check interpolation to TRANSP grid
    figure;
    plot(sqrt_psin_chain2,EL.NE(ijp,:)-6.0*ION.NC(ijp,:),'bx')
    hold on
    plot(sqrt_psin_TRANSP,ni,'r-')    
    legend('JETPEAK chain2','fit to TRANSP')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$n_i$','Interpreter','LaTex')
    % Then check derivative
    figure;
    plot(rpsi_TRANSP,ni,'b-x')
    hold on
    plot(rpsi_TRANSP,...
        ni(iflxsurf)+...
           (rpsi_TRANSP-rpsi_TRANSP(iflxsurf))*dni_drho(iflxsurf),'r-')
    legend('Experimental data','Tangent at point of interest')
    xlabel('$r_\psi$','Interpreter','LaTex')
    ylabel('$n_i$','Interpreter','LaTex')
    hold on
    plot(rpsi_TRANSP(iflxsurf),ni(iflxsurf),'ro')
    hold off
    
    % Electron density
    % First check interpolation to TRANSP grid
    figure;
    plot(sqrt_psin_chain2,EL.NE(ijp,:),'bx')
    hold on
    plot(sqrt_psin_TRANSP,ne,'r-')    
    legend('JETPEAK chain2','GS2 fit')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$n_e$','Interpreter','LaTex')
    % Then check derivative
    figure;
    plot(rpsi_TRANSP,ne,'b-x')
    hold on
    plot(rpsi_TRANSP, ...
        ne(iflxsurf)+...
           (rpsi_TRANSP-rpsi_TRANSP(iflxsurf))*dne_drho(iflxsurf),'r-')
    legend('Experimental data','Tangent at point of interest')
    xlabel('$r_\psi$','Interpreter','LaTex')
    ylabel('$n_e$','Interpreter','LaTex')
    hold on
    plot(rpsi_TRANSP(iflxsurf),ne(iflxsurf),'ro')
    hold off
    
    % Carbon density
    % First check interpolation to TRANSP grid
    figure;
    plot(sqrt_psin_chain2,ION.NC(ijp,:),'bx')
    hold on
    plot(sqrt_psin_TRANSP,nc,'r-')    
    legend('JETPEAK chain2','GS2 fit')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$n_C$','Interpreter','LaTex')
    % Then check derivative
    figure;
    plot(rpsi_TRANSP,nc,'b-x')
    hold on
    plot(rpsi_TRANSP,...
        nc(iflxsurf)+...
           (rpsi_TRANSP-rpsi_TRANSP(iflxsurf))*dnc_drho(iflxsurf),'r-')
    legend('Experimental data','Tangent at point of interest')
    xlabel('$r_\psi$','Interpreter','LaTex')
    ylabel('$n_C$','Interpreter','LaTex')
    hold on
    plot(rpsi_TRANSP(iflxsurf),nc(iflxsurf),'ro')
    hold off
    
    % Ion temperature
    % First check interpolation to TRANSP grid
    figure;
    plot(sqrt_psin_chain2,ION.TI(ijp,:),'bx')
    hold on
    plot(sqrt_psin_TRANSP,ti,'r-')    
    legend('JETPEAK chain2','GS2 fit')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$T_i$','Interpreter','LaTex')
    % Then check derivative
    figure;
    plot(rpsi_TRANSP,ti,'b-x')
    hold on
    plot(rpsi_TRANSP,...
        ti(iflxsurf)+...
           (rpsi_TRANSP-rpsi_TRANSP(iflxsurf))*dti_drho(iflxsurf),'r-')
    legend('Experimental data','Tangent at point of interest')
    xlabel('$r_\psi$','Interpreter','LaTex')
    ylabel('$T_i$','Interpreter','LaTex')
    hold on
    plot(rpsi_TRANSP(iflxsurf),ti(iflxsurf),'ro')
    hold off
    
    % Electron temperature
    % First check interpolation to TRANSP grid
    figure;
    plot(sqrt_psin_chain2,EL.TE(ijp,:),'bx')
    hold on
    plot(sqrt_psin_TRANSP,te,'r-')    
    legend('JETPEAK chain2','GS2 fit')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$T_e$','Interpreter','LaTex')
    % Then check derivative
    figure;
    plot(rpsi_TRANSP,te,'b-x')
    hold on
    plot(rpsi_TRANSP,...
        te(iflxsurf)+...
           (rpsi_TRANSP-rpsi_TRANSP(iflxsurf))*dte_drho(iflxsurf),'r-')
    legend('Experimental data','Tangent at point of interest')
    xlabel('$r_\psi$','Interpreter','LaTex')
    ylabel('$T_e$','Interpreter','LaTex')
    hold on
    plot(rpsi_TRANSP(iflxsurf),te(iflxsurf),'ro')
    hold off
    
    % Carbon temperature: redundant, since T_C = T_i above.
    
    % Plasma pressure
    figure;
    plot(rpsi_TRANSP,p,'b-x')
    hold on
    plot(rpsi_TRANSP,...
        p(iflxsurf)+...
           (rpsi_TRANSP-rpsi_TRANSP(iflxsurf))*dp_drho(iflxsurf),'r-')
    legend('Experimental data','Tangent at point of interest')
    xlabel('$r_\psi$','Interpreter','LaTex')
    ylabel('$p$','Interpreter','LaTex')
    hold on
    plot(rpsi_TRANSP(iflxsurf),p(iflxsurf),'ro')
    hold off
    
    % Collisionality
    figure;
    plot(sqrt_psin_TRANSP,loglamda,'r-')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$\lambda_{ei}$','Interpreter','LaTex')
    figure;
    plot(sqrt_psin_TRANSP,nu_ii,'r-')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$\nu_{ii}$','Interpreter','LaTex')
    figure;
    plot(sqrt_psin_TRANSP,nu_ee,'r-')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$\nu_{ee}$','Interpreter','LaTex')
    
    % Toroidal angular velocity
    % First check interpolation to TRANSP grid
    figure;
    plot(sqrt_psin_chain2,abs(ION.ANGF(ijp,:)),'bx')
    hold on
    plot(sqrt_psin_TRANSP,abs(omega),'r-')    
    legend('JETPEAK chain2','GS2 fit')
    xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
    ylabel('$\vert\Omega_\zeta\vert$','Interpreter','LaTex')
    title(['sgn(Omega)= ' num2str(sign_mach)], 'FontSize',18)
    % Then check derivative
    figure;
    plot(rpsi_TRANSP,omega,'b-x')
    hold on
    plot(rpsi_TRANSP,...
        omega(iflxsurf)+...
           (rpsi_TRANSP-rpsi_TRANSP(iflxsurf))*domega_drho(iflxsurf),'r-')
    legend('Experimental data','Tangent at point of interest')
    xlabel('$r_\psi$','Interpreter','LaTex')
    ylabel('$\Omega_\zeta$','Interpreter','LaTex')
    hold on
    plot(rpsi_TRANSP(iflxsurf),omega(iflxsurf),'ro')
    hold off
    
end


end

