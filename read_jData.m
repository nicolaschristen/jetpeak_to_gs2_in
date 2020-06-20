%% Read plasma parameters and geometry from JETPEAK, TRANSP and ASCOT databases
%
% Input :   ijp -- shot index in JETPEAK DB
%           skip_I -- [kw, 1]when set to 1, skips computation of I(psi)
%           check_iFlxSurf -- [kw, []] list of flux srufaces to plot if visual
%                             check is needed
%           showWarnings -- [kw, 0] show warnings when part of the data is
%                           incomplete
%           trinity_norm -- [kw, 0] if true, gs2 flux dotted with gradPsi
%                           else dotted with grad(x).
%
% Output:   jData -- structure with dimensionful extracted data
%
% All dimensionful quantities are in SI units,
% except for temperatures which are in eV.
%
function jData = read_jData(ijp, varargin)

opt_defaults = struct( 'skip_I', 0, ...
                       'check_iFlxSurf', [], ...
                       'showWarnings', 0, ...
                       'trinity_norm', 0 );
opt = get_optargin(opt_defaults, varargin);

% loading databases
load ~/codes/jetpeak_v_gs2/databases/TRANSP_2017_3.mat
load ~/codes/jetpeak_v_gs2/databases/JETPEAK_2019_10.mat

% Ignore NaN warnings
warning('off','MATLAB:chckxy:IgnoreNaN');

% Indices corresponding to ijp
itransp=find(TRANSP.JPI==ijp);
idxQAscot = find(Q.ASCOT.JPI==ijp);
idxAscot = find(ASCOT.SAMPLE==ijp);
if ~isempty(idxAscot)
    % Per ijp, there are many ASCOT runs with slightly different
    % settings. Here, we simply pick the first run in the list.
    idxAscot = idxAscot(1);

% JET shot number corresponding to ijp
shot = SAMPLE.SHOT(ijp);

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
sign_PI_ASCOT_v_GS2 = -1;




%% Extract mag. geometry from TRANSP

% Poloidal flux [T m^2]
% with definition equivalent to GS2's, i.e. without toroidal 2pi factor.
% See Edmund's thesis and the TRANSP variable list.
psi=permute(TRANSP.T.PSI,[3 2 1]); % on rectangular grid
psi=psi(:,:,itransp); % psi(iR,iZ)
psiflu=permute(TRANSP.T.PLFLX,[2 1]); % on flux surface grid
psiflu=psiflu(:,itransp); 
psiflu=psiflu.';
nflxsurf=numel(psiflu); % number of flux surf.

% Volume elements dV(i) = V(psi(i))-V(psi(i-1))
dV = TRANSP.G.DVOL(itransp,:);

% Volume of flux surfaces
V = zeros(1, nflxsurf);
V(1) = dV(1);
for iFlxSurf = 2:nflxsurf
    V(iFlxSurf) = V(iFlxSurf-1) + dV(iFlxSurf);
end

% dV/dpsi
derivOrder = 1;
dV_dpsi = interpol(psiflu, V, psiflu, derivOrder);
    
% Flux surface areas from TRANSP [m^{-2}]
A_psi = 1e-4*TRANSP.T.SURF(itransp,:);

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
rpsi_TRANSP=zeros(1,nflxsurf); % GS2 definition of rpsi for irho=2, not yet normalized
Rmaj=zeros(1,nflxsurf); % Rmaj definition for iflux ~= 1, not yet normalized
Rmax = TRANSP.G.RMAJM(itransp,2*nflxsurf+1); % Maximum major radius of confined plasma region
for indx=1:nflxsurf
    rpsi_TRANSP(indx)= ...
        (TRANSP.G.RMAJM(itransp,nflxsurf+1+indx) ...
           - TRANSP.G.RMAJM(itransp,nflxsurf+1-indx))/2.;
    Rmaj(indx)=TRANSP.G.RMAJM(itransp,nflxsurf+1+indx)-rpsi_TRANSP(indx);
end
% Free choice for Rgeo. This then defines Bref., see notes about signs
Rgeo=sign_Rgeo*abs(Rmaj);
dR_drpsi=interpol(rpsi_TRANSP,Rmaj,rpsi_TRANSP,1);

% Vertical coord. [m]
Zmag=1.e-2*TRANSP.T.YAXIS(itransp);
Zflu=1.e-2*permute(TRANSP.T.ZFLU,[2 3 1]);
Zflu=Zflu(:,:,itransp); % (iflx,itheta)
Zpsi= 1.e-2*TRANSP.T.PSIZ;

% Compute <|grad(psi)|>_psi
gradPsiAvg = compute_gradPsiAvg(Rpsi,Zpsi,psi,Rflu,Zflu,Rmag,Zmag);

% Elongation []
kappa=permute(TRANSP.T.ELONG,[2 1]);
kappa=kappa(:,itransp);
kappa=kappa.';
dkappa_drpsi=interpol(rpsi_TRANSP,kappa,rpsi_TRANSP,1);

% Triangularity []
delta=permute(TRANSP.T.TRIANG,[2 1]);
delta=asin(delta(:,itransp));
delta=delta.';
ddelta_drpsi=interpol(rpsi_TRANSP,delta,rpsi_TRANSP,1);




%% Extract q and plasma profiles from EFIT, EL and ION

% Normalized sqrt(psi) grid used in chain2 []
sqrt_psin_chain2=linspace(0.,1.,21);

% Safety factor []
q = sign_q * interpol(EFIT.RMJO(ijp,:), ...
    abs(EFIT.Q(ijp,:)),TRANSP.G.RMAJM(itransp,nflxsurf+2:end));
dq_drpsi=interpol(rpsi_TRANSP,q,rpsi_TRANSP,1);
shat = rpsi_TRANSP./q .* dq_drpsi;

% Density [m^{-3}]
ni=EL.NE(ijp,:)-6.0*ION.NC(ijp,:);
ni=interpol(sqrt_psin_chain2,ni,sqrt_psin_TRANSP);
nc=interpol(sqrt_psin_chain2,ION.NC(ijp,:),sqrt_psin_TRANSP);
ne=interpol(sqrt_psin_chain2,EL.NE(ijp,:),sqrt_psin_TRANSP);
% Should carbon be included?
add_carbon = zeros(1,nflxsurf);
for iFlxSurf = 1:nflxsurf
    if nc(iFlxSurf)/ni(iFlxSurf) >= 0.05
        % yes
        add_carbon(iFlxSurf) = 1;
    else
        % if not, then adapt ni so that QN holds
        add_carbon(iFlxSurf) = 0;
        ni(iFlxSurf) = ni(iFlxSurf) + Zimp/Zmain*nc(iFlxSurf);
    end
end
dni_drpsi=interpol(rpsi_TRANSP,ni,rpsi_TRANSP,1);
dne_drpsi=interpol(rpsi_TRANSP,ne,rpsi_TRANSP,1);
dnc_drpsi=interpol(rpsi_TRANSP,nc,rpsi_TRANSP,1);

% Temperature [eV]
ti=interpol(sqrt_psin_chain2,ION.TI(ijp,:),sqrt_psin_TRANSP);
dti_drpsi=interpol(rpsi_TRANSP,ti,rpsi_TRANSP,1);
te=interpol(sqrt_psin_chain2,EL.TE(ijp,:),sqrt_psin_TRANSP);
dte_drpsi=interpol(rpsi_TRANSP,te,rpsi_TRANSP,1);
tc=ti; % no direct measurement for T_imp
dtc_drpsi=dti_drpsi;

% Plasma pressure [Pa]
p=e*ni.*ti+e*ne.*te;
p(add_carbon==1) = p(add_carbon==1)+e*nc(add_carbon==1).*tc(add_carbon==1);
dp_drpsi=interpol(rpsi_TRANSP,p,rpsi_TRANSP,1);

% Collisionality (See notes) [s^{-1}]
colfac = e^(-3/2)/(4*pi*eps0)^2; % conversion factor from cgs to SI+eV
Zi = 1; % main ion is deuterium
Ze = -1;
mi = 2*mp; % main ion is deuterium
loglamda = coulomb_log(Zi,mi,ne,ni,te,ti);
nu_ii = colfac * 4*pi*Zi^4*e^4*ni.*loglamda./(sqrt(mi)*(2*ti).^(3/2));
nu_ee = colfac * 4*pi*Ze^4*e^4*ne.*loglamda./(sqrt(me)*(2*te).^(3/2));
Zc = 6;
mc = 12*mp;
loglamda_c = coulomb_log(Zc,mc,ne,nc,te,tc);
nu_cc = colfac * 4*pi*Zc^4*e^4*nc.*loglamda_c./(sqrt(mc)*(2*tc).^(3/2));

% Tor. angular vel. [radian/s]
omega = interpol(sqrt_psin_chain2,ION.ANGF(ijp,:),sqrt_psin_TRANSP);
domega_drpsi = interpol(rpsi_TRANSP,omega,rpsi_TRANSP,1);
% see notes about signs
sign_mach = sign(omega)*sign(BASIC.IP(ijp));
% see notes about signs
sign_domega_drpsi = sign(domega_drpsi)*sign(BASIC.IP(ijp));
omega = sign_mach.*abs(omega);
domega_drpsi = sign_domega_drpsi.*abs(domega_drpsi);

% Effective charge [C]
Zeff = zeros(1,nflxsurf);
for iFlxSurf = 1:nflxsurf
    if add_carbon(iFlxSurf)
       Zeff(iFlxSurf) = (Zmain^2*ni(iFlxSurf)+Zimp^2*nc(iFlxSurf))./(Zmain*ni(iFlxSurf)+Zimp*nc(iFlxSurf));
    else
       Zeff(iFlxSurf) = Zmain;
    end
end




%% Calculate I(psi) to determine Bref for given Rgeo
% where I(psi) comes from expressing B as B = I*grad(zeta) + grad(zeta) x grad(psi)

I = zeros(1,nflxsurf);

if ~opt.skip_I

    for iFlxSurf = 1:nflxsurf
        
        if ismember(iFlxSurf, opt.check_iFlxSurf)
            plot_verbose_Ipsi_integration_loc = plot_verbose_Ipsi_integration;
        else
            plot_verbose_Ipsi_integration_loc = 0;
        end
        
        Rflu_loc=permute(Rflu,[2 1]);
        Rflu_loc=Rflu_loc(:,iFlxSurf); % Rflu for all poloidal locations on chosen flux surf
        Zflu_loc=permute(Zflu,[2 1]);
        Zflu_loc=Zflu_loc(:,iFlxSurf); % Zflu for all poloidal locations on chosen flux surf
        
        try
            I(iFlxSurf) = compute_I(Rmag,Zmag,Rflu_loc,Zflu_loc,Rpsi,Zpsi,psi,q(iFlxSurf), ...
                psiflu(iFlxSurf), plot_verbose_Ipsi_integration_loc);
        catch err
            I(iFlxSurf) = NaN;
            if opt.showWarnings
                fprintf('\nCannot compute I(psi) for flux surface index %d\nError message: %s\nSkipping it.\n\n',iFlxSurf,err.identifier)
            end
        end
        
    end % loop over flux surfaces

end

% In GS2, Bref is defined via Bref=I/(Rgeo*a)
Bref=abs(I./Rgeo); % no 1/a since Rgeo not normalised yet, see notes about signs




%% Compute experimental heat fluxes from PENCIL and ASCOT

% Conversion factor to GS2 radial coordinate
dx_dpsi = abs(q./(rpsi_TRANSP.*Bref));

% Energy deposition from PENCIL [W/m^3]
srcE_i_PENCIL = interpol(sqrt_psin_chain2, Q.NBPICH.qis(ijp,:), ...
    sqrt_psin_TRANSP); % ions
srcE_e_PENCIL = interpol(sqrt_psin_chain2, Q.NBPICH.qes(ijp,:), ...
    sqrt_psin_TRANSP); % e-
srcE_ie_PENCIL = interpol(sqrt_psin_chain2, Q.NBPICH.qie(ijp,:), ...
    sqrt_psin_TRANSP); % ions to e- transfer
% Incorporate transfer term
srcE_i_PENCIL = srcE_i_PENCIL - srcE_ie_PENCIL;
srcE_e_PENCIL = srcE_e_PENCIL + srcE_ie_PENCIL;
% Energy deposition from Q.ASCOT [W/m^3]
if isempty(idxQAscot)
    fprintf('JETPEAK index %d has no corresponding index in Q.ASCOT.',ijp)
    srcE_i_QASCOT = NaN(nflxsurf,1);
    srcE_e_QASCOT = NaN(nflxsurf,1);
    srcE_ie_QASCOT = NaN(nflxsurf,1);
else
    srcE_i_QASCOT = interpol(sqrt_psin_chain2, Q.ASCOT.qis(idxQAscot,:), ...
        sqrt_psin_TRANSP); % ions
    srcE_e_QASCOT = interpol(sqrt_psin_chain2, Q.ASCOT.qes(idxQAscot,:), ...
        sqrt_psin_TRANSP); % e-
    srcE_ie_QASCOT = interpol(sqrt_psin_chain2, Q.ASCOT.qie(idxQAscot,:), ...
        sqrt_psin_TRANSP); % ions to e- transfer
end
% Incorporate transfer term
srcE_i_QASCOT = srcE_i_QASCOT - srcE_ie_QASCOT;
srcE_e_QASCOT = srcE_e_QASCOT + srcE_ie_QASCOT;

% Associated ion and electron heat fluxes (see notes) [kg/(s^3)]
Qi_PENCIL = flux_from_source(psiflu, dV, srcE_i_PENCIL);
Qe_PENCIL = flux_from_source(psiflu, dV, srcE_e_PENCIL);
Qi_QASCOT = flux_from_source(psiflu, dV, srcE_i_QASCOT);
Qe_QASCOT = flux_from_source(psiflu, dV, srcE_e_QASCOT);




%% Compute experimental toroidal angular momentum flux from ASCOT and PENCIL

% Torque deposition from PENCIL [Pa]
if isnan(NB.NBP.NBP2TORP(ijp,1))
    if opt.showWarnings
        fprintf('\nJETPEAK index %d has no data in NB.NBP.NBP2TORP.\n',ijp)
    end
    srcL_PENCIL = NaN(nflxsurf,1);
else
    srcL_PENCIL = interpol(sqrt_psin_chain2, NB.NBP.NBP2TORP(ijp,:), sqrt_psin_TRANSP);
end
% Torque deposition from ASCOT [Pa]
if isempty(idxAscot)
    fprintf('\nJETPEAK index %d has no corresponding index in ASCOT.\n',ijp)
    srcL_ASCOT = NaN(nflxsurf,1);
else
    srcL_ASCOT = interpol( ...
        sqrt_psin_chain2, ...
        ASCOT.COLLTQI(idxAscot,:) + ASCOT.COLLTQE(idxAscot,:) ...
            + ASCOT.COLLTQIMP(idxAscot,:) + ASCOT.JXBTORQ(idxAscot,:), ...
        sqrt_psin_TRANSP);
    % Sign conventions in GS2 and ASCOT are opposite:
    % in GS2, negative = inward,
    % in ASCOT, negative = outward.
    % Take this sign into account here:
    srcL_ASCOT = sign_PI_ASCOT_v_GS2 * srcL_ASCOT;
end

% Associated toroidal angular momentum flux (see notes) [kg/s^2]
PI_PENCIL = flux_from_source(psiflu, dV, dx_dpsi, srcL_PENCIL);
PI_ASCOT = flux_from_source(psiflu, dV, dx_dpsi, srcL_ASCOT);




%% Collect into jData structure

jData.a=a;
jData.nref=ni;
jData.tref=ti;
jData.mref=1.675e-27+1.6726e-27; % assume main ion is deuterium
jData.Bref=Bref;
jData.vthref=sqrt(2.0*e*jData.tref/jData.mref);
jData.rhoref=sqrt(2.0*jData.mref*e*jData.tref)./(e*Bref);
jData.rhostar=jData.rhoref/a;

jData.shot = shot;
jData.ijp = ijp;
jData.idxQAscot = idxQAscot;
jData.idxAscot = idxAscot;
jData.Rmag=Rmag;
jData.Zmag=Zmag;
jData.Rmaj=Rmaj;
jData.Rmax=Rmax;
jData.Rgeo=Rgeo;
jData.sqrt_psin_TRANSP=sqrt_psin_TRANSP;
jData.sqrt_psin_chain2=sqrt_psin_chain2;
jData.dV=dV;
jData.V=V;
jData.dV_dpsi=dV_dpsi;
jData.A_psi=A_psi;
jData.rpsi=rpsi_TRANSP;
jData.psiflu = psiflu;
jData.dx_dpsi = dx_dpsi;
jData.dR_drpsi=dR_drpsi;
jData.gradPsiAvg = gradPsiAvg;
jData.q=q;
jData.shat=shat;
jData.kappa=kappa;
jData.dkappa_drpsi=dkappa_drpsi;
jData.delta=delta;
jData.ddelta_drpsi=ddelta_drpsi;

jData.ni=ni;
jData.dni_drpsi=dni_drpsi;
jData.ne=ne;
jData.dne_drpsi=dne_drpsi;
jData.add_carbon=add_carbon;
jData.nc=nc;
jData.dnc_drpsi=dnc_drpsi;
jData.ti=ti;
jData.dti_drpsi=dti_drpsi;
jData.te=te;
jData.dte_drpsi=dte_drpsi;
jData.tc=tc;
jData.dtc_drpsi=dtc_drpsi;
jData.p=p;
jData.dp_drpsi=dp_drpsi;
jData.nu_ii=nu_ii;
jData.nu_ee=nu_ee;
jData.nu_cc=nu_cc;
jData.omega=omega;
jData.domega_drpsi=domega_drpsi;
jData.Zeff=Zeff;

jData.srcE_i_PENCIL = srcE_i_PENCIL;
jData.srcE_i_QASCOT = srcE_i_QASCOT;
jData.srcE_e_PENCIL = srcE_e_PENCIL;
jData.srcE_e_QASCOT = srcE_e_QASCOT;
jData.srcE_ie_PENCIL = srcE_ie_PENCIL;
jData.srcE_ie_QASCOT = srcE_ie_QASCOT;
jData.Qi_PENCIL = Qi_PENCIL;
jData.Qi_QASCOT = Qi_QASCOT;
jData.Qe_PENCIL = Qe_PENCIL;
jData.Qe_QASCOT = Qe_QASCOT;
% Normalisation factor
if opt.trinity_norm
    jData.QNorm = jData.nref*e.*jData.tref.*jData.vthref.*jData.rhostar.^2.*jData.gradPsiAvg;
else
    jData.QNorm = jData.nref*e.*jData.tref.*jData.vthref.*jData.rhostar.^2./jData.dx_dpsi;
end

jData.srcL_PENCIL = srcL_PENCIL;
jData.srcL_ASCOT = srcL_ASCOT;
jData.PI_PENCIL = PI_PENCIL;
jData.PI_ASCOT = PI_ASCOT;
% Normalisation factor
if opt.trinity_norm
    jData.PINorm = jData.nref*jData.mref.*jData.vthref.^2.*jData.a.*jData.rhostar.^2 ...
                   .*jData.gradPsiAvg;
else
    jData.PINorm = jData.nref*jData.mref.*jData.vthref.^2.*jData.a.*jData.rhostar.^2 ...
                   ./jData.dx_dpsi;
end

% Normalisation factor for the particle flux
if opt.trinity_norm
    jData.GammaNorm = jData.nref.*jData.vthref.*jData.rhostar.^2.*jData.gradPsiAvg;
else
    jData.GammaNorm = jData.nref.*jData.vthref.*jData.rhostar.^2./jData.dx_dpsi;
end




%% Check data prof., interpolations, derivatives (if plot_verbose)

for iplot = 1:numel(opt.check_iFlxSurf)

    iFlxSurf = opt.check_iFlxSurf(iplot);
    
    % Major radius
    figure;
    plot(rpsi_TRANSP,Rmaj,'b-x')
    hold on
    plot(rpsi_TRANSP, ...
        Rmaj(iFlxSurf)+ ...
           (rpsi_TRANSP-rpsi_TRANSP(iFlxSurf))*dR_drpsi(iFlxSurf),'r-')
    legend('Experimental data','Tangent at point of interest')
    xlabel('$r_\psi$','Interpreter','LaTex')
    ylabel('$R_\psi$','Interpreter','LaTex')
    hold on
    plot(rpsi_TRANSP(iFlxSurf),Rmaj(iFlxSurf),'ro')
    hold off
    
    % Elongation
    figure;
    plot(rpsi_TRANSP,kappa,'b-x')
    hold on
    plot(rpsi_TRANSP, ...
        kappa(iFlxSurf)+ ...
           (rpsi_TRANSP-rpsi_TRANSP(iFlxSurf))*dkappa_drpsi(iFlxSurf),'r-')
    legend('Experimental data','Tangent at point of interest')
    xlabel('$r_\psi$','Interpreter','LaTex')
    ylabel('$\kappa$','Interpreter','LaTex')
    hold on
    plot(rpsi_TRANSP(iFlxSurf),kappa(iFlxSurf),'ro')
    hold off
    
    % Triangularity
    figure;
    plot(rpsi_TRANSP,delta,'b-x')
    hold on
    plot(rpsi_TRANSP, ...
        delta(iFlxSurf)+...
           (rpsi_TRANSP-rpsi_TRANSP(iFlxSurf))*ddelta_drpsi(iFlxSurf),'r-')
    legend('Experimental data','Tangent at point of interest')
    xlabel('$r_\psi$','Interpreter','LaTex')
    ylabel('$\delta$','Interpreter','LaTex')
    hold on
    plot(rpsi_TRANSP(iFlxSurf),delta(iFlxSurf),'ro')
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
        q(iFlxSurf)+...
           (rpsi_TRANSP-rpsi_TRANSP(iFlxSurf))*dq_drpsi(iFlxSurf),'r-')
    legend('Experimental data','Tangent at point of interest')
    xlabel('$r_\psi$','Interpreter','LaTex')
    ylabel('$q$','Interpreter','LaTex')
    hold on
    plot(rpsi_TRANSP(iFlxSurf),q(iFlxSurf),'ro')
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
        ni(iFlxSurf)+...
           (rpsi_TRANSP-rpsi_TRANSP(iFlxSurf))*dni_drpsi(iFlxSurf),'r-')
    legend('Experimental data','Tangent at point of interest')
    xlabel('$r_\psi$','Interpreter','LaTex')
    ylabel('$n_i$','Interpreter','LaTex')
    hold on
    plot(rpsi_TRANSP(iFlxSurf),ni(iFlxSurf),'ro')
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
        ne(iFlxSurf)+...
           (rpsi_TRANSP-rpsi_TRANSP(iFlxSurf))*dne_drpsi(iFlxSurf),'r-')
    legend('Experimental data','Tangent at point of interest')
    xlabel('$r_\psi$','Interpreter','LaTex')
    ylabel('$n_e$','Interpreter','LaTex')
    hold on
    plot(rpsi_TRANSP(iFlxSurf),ne(iFlxSurf),'ro')
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
        nc(iFlxSurf)+...
           (rpsi_TRANSP-rpsi_TRANSP(iFlxSurf))*dnc_drpsi(iFlxSurf),'r-')
    legend('Experimental data','Tangent at point of interest')
    xlabel('$r_\psi$','Interpreter','LaTex')
    ylabel('$n_C$','Interpreter','LaTex')
    hold on
    plot(rpsi_TRANSP(iFlxSurf),nc(iFlxSurf),'ro')
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
        ti(iFlxSurf)+...
           (rpsi_TRANSP-rpsi_TRANSP(iFlxSurf))*dti_drpsi(iFlxSurf),'r-')
    legend('Experimental data','Tangent at point of interest')
    xlabel('$r_\psi$','Interpreter','LaTex')
    ylabel('$T_i$','Interpreter','LaTex')
    hold on
    plot(rpsi_TRANSP(iFlxSurf),ti(iFlxSurf),'ro')
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
        te(iFlxSurf)+...
           (rpsi_TRANSP-rpsi_TRANSP(iFlxSurf))*dte_drpsi(iFlxSurf),'r-')
    legend('Experimental data','Tangent at point of interest')
    xlabel('$r_\psi$','Interpreter','LaTex')
    ylabel('$T_e$','Interpreter','LaTex')
    hold on
    plot(rpsi_TRANSP(iFlxSurf),te(iFlxSurf),'ro')
    hold off
    
    % Carbon temperature: redundant, since T_C = T_i above.
    
    % Plasma pressure
    figure;
    plot(rpsi_TRANSP,p,'b-x')
    hold on
    plot(rpsi_TRANSP,...
        p(iFlxSurf)+...
           (rpsi_TRANSP-rpsi_TRANSP(iFlxSurf))*dp_drpsi(iFlxSurf),'r-')
    legend('Experimental data','Tangent at point of interest')
    xlabel('$r_\psi$','Interpreter','LaTex')
    ylabel('$p$','Interpreter','LaTex')
    hold on
    plot(rpsi_TRANSP(iFlxSurf),p(iFlxSurf),'ro')
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
        omega(iFlxSurf)+...
           (rpsi_TRANSP-rpsi_TRANSP(iFlxSurf))*domega_drpsi(iFlxSurf),'r-')
    legend('Experimental data','Tangent at point of interest')
    xlabel('$r_\psi$','Interpreter','LaTex')
    ylabel('$\Omega_\zeta$','Interpreter','LaTex')
    hold on
    plot(rpsi_TRANSP(iFlxSurf),omega(iFlxSurf),'ro')
    hold off
    
end


end

