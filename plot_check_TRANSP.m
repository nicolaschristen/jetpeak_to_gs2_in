%% Plot quantities to compare between different diagnostics
%
% Input : itransp -- shot index in TRANSP
%         ijp -- shot index in JETPEAK
%         EFIT -- structure obtained by loading JETPEAK DB
%         EL -- structure obtained by loading JETPEAK DB
%         ION -- structure obtained by loading JETPEAK DB
%         TRANSP -- structure obtained by loading TRANSP DB
%
% Output: -
%
% According to experimentalists, JETPEAK data is
% more reliable than TRANSP data.
%
function plot_check_TRANSP(itransp,ijp,EFIT,EL,ION,TRANSP)

% Need access to ../interpol.m
addpath(genpath('..'))

% chain 2 radial coordinate
rho_chain2=linspace(0.,1.,21);

% get normalized sqrt(psi) for TRANSP data
nflxsurf=size(TRANSP.T.PLFLX,2);
rho_TRANSP=zeros(1,nflxsurf);
for i=1:nflxsurf
    rho_TRANSP(i)=sqrt(TRANSP.T.PLFLX(itransp,i)/TRANSP.T.PLFLX(itransp,end));
end

% get psi(R,Z) matrix for EFIT data
EPSI=zeros(numel(EFIT.PSIR),numel(EFIT.PSIZ));
for i=1:numel(EFIT.PSIZ)
    for j=0:numel(EFIT.PSIR)-1
        EPSI(j+1,i)=EFIT.PSI(ijp,33*j+i);
    end
end

% get R,Z from TRANSP
TR_flu=1.e-2*permute(TRANSP.T.RFLU,[2 3 1]); % flux surface grid
TR_flu=TR_flu(:,:,itransp);
TR_rect= 1.e-2*TRANSP.T.PSIR; % rectangular grid
TZ_flu=1.e-2*permute(TRANSP.T.ZFLU,[2 3 1]);
TZ_flu=TZ_flu(:,:,itransp);
TZ_rect= 1.e-2*TRANSP.T.PSIZ;
ntheta = size(TR_flu,2);

% get psi(R,Z) matrix for TRANSP data
TPSI_rect=permute(TRANSP.T.PSI,[3 2 1]); % on rectangular grid
TPSI_rect=TPSI_rect(:,:,itransp); % psi(iR,iZ)
TPSI_flu=permute(TRANSP.T.PLFLX,[2 1]); % on flux surface grid
TPSI_flu=TPSI_flu(:,itransp); % psiflu(iflxsurf)
TPSI_flu_thet = zeros(nflxsurf,ntheta);
for iflx=1:nflxsurf
    TPSI_flu_thet(iflx,:) = TPSI_flu(iflx);
end
% Toroidal mag flux
TPSItor_flu=permute(TRANSP.T.TRFLX,[2 1]); % on flux surface grid
TPSItor_flu=TPSItor_flu(:,itransp); % psiflu(iflxsurf)
TPSItor_flu_thet = zeros(nflxsurf,ntheta);
for iflx=1:nflxsurf
    TPSItor_flu_thet(iflx,:) = TPSItor_flu(iflx);
end

% interpolate q-profile from EFIT to TRANSP radial coord.
EQ=interpol(EFIT.RMJO(ijp,:),EFIT.Q(ijp,:),TRANSP.G.RMAJM(itransp,nflxsurf+2:end));

% Compare Te
figure
plot(rho_chain2,EL.TE(ijp,:),'b-x')
hold on
plot(rho_TRANSP,TRANSP.G.TE(itransp,:),'r-*')
yl=ylim(gca);
ylim(gca,[0,yl(2)]);
xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
ylabel('$T_e$ [eV]','Interpreter','LaTex')
legend('HRTS+LIDAR','TRANSP')
% Compare ne
figure
plot(rho_chain2,EL.NE(ijp,:),'b-x')
hold on
plot(rho_TRANSP,TRANSP.G.NE(itransp,:),'r-*')
yl=ylim(gca);
ylim(gca,[0,yl(2)]);
xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
ylabel('$n_e$ [m**-3]','Interpreter','LaTex')
legend('HRTS+LIDAR','TRANSP')
% Compare Ti
figure
plot(rho_chain2,ION.TI(ijp,:),'b-x')
hold on
plot(rho_TRANSP,TRANSP.G.TI(itransp,:),'r-*')
yl=ylim(gca);
ylim(gca,[0,yl(2)]);
xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
ylabel('$T_i$ [eV]','Interpreter','LaTex')
legend('CXS','TRANSP')
% Compare ni
figure
plot(rho_chain2,EL.NE(ijp,:)-6.*ION.NC(ijp,:),'b-x')
hold on
plot(rho_TRANSP,TRANSP.G.ND(itransp,:),'r-*')
yl=ylim(gca);
ylim(gca,[0,yl(2)]);
xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
ylabel('$n_i$ [m**-3]','Interpreter','LaTex')
legend('HRTS+LIDAR','TRANSP')
% Compare omega
figure
plot(rho_chain2,ION.ANGF(ijp,:),'b-x')
hold on
plot(rho_TRANSP,TRANSP.G.OMEGA(itransp,:),'r-*')
yl=ylim(gca);
ylim(gca,[0,yl(2)]);
xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
ylabel('$\omega$ [s**-1]','Interpreter','LaTex')
legend('CXS','TRANSP')
% Compare q
figure
plot(rho_TRANSP,EQ,'b-x')
hold on
plot(rho_TRANSP,TRANSP.T.Q(itransp,:),'r-*')
yl=ylim(gca);
ylim(gca,[0,yl(2)]);
xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
ylabel('$q$','Interpreter','LaTex')
legend('EFIT','TRANSP')
% Compare flux surface shape
figure
contourf(TR_flu,TZ_flu,TPSI_flu_thet)
hold on
contour(EFIT.PSIR,EFIT.PSIZ,EPSI,'LineColor','r','LineWidth',2)
colorbar
caxis([min(TPSI_flu) max(TPSI_flu)])
title('$\psi$ from TRANSP flux-surface grid, red is EFIT')
figure
contourf(TR_rect,TZ_rect,TPSI_rect')
hold on
contour(EFIT.PSIR,EFIT.PSIZ,EPSI,'LineColor','r','LineWidth',2)
colorbar
caxis([min(min(TPSI_rect)) max(max(TPSI_rect))])
title('$\psi$ from TRANSP rectangular grid, red is EFIT')
% Compare toroidal magnetic flux
figure('Position',[10,10,1500,600])
subplot(1,2,1)
contourf(TR_flu,TZ_flu,TPSI_flu_thet)
myxlim = xlim();
myylim = ylim();
colorbar
mycaxis = caxis;
title('Tor mag flux from TRANSP')
subplot(1,2,2)
contourf(EFIT.PSIR,EFIT.PSIZ,EPSI/6)
xlim(myxlim)
ylim(myylim)
colorbar
caxis(mycaxis)
title('Tor mag flux from EFIT')

% Plot gradient length of Te, should get : rise-plateau-rise-fall
dTe_drho=interpol(rho_chain2,EL.TE(ijp,:),rho_chain2,1);
figure
plot(rho_chain2,-1.*dTe_drho./EL.TE(ijp,:),'b-x')
xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
ylabel('$\vert \nabla(\ln(T_e)) \vert$ [m**-1]','Interpreter','LaTex')
legend('HRTS+LIDAR')


end
