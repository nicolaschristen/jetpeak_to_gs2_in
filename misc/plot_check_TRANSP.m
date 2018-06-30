function plot_check_TRANSP(itransp,ijp,EFIT,EL,ION,TRANSP)

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

% get psi(R,Z) matrix for TRANSP data
TPSI=permute(TRANSP.T.PSI,[2 3 1]);
TPSI=TPSI(:,:,itransp);

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
% Compare psi
figure
contour(EFIT.PSIR,EFIT.PSIZ,EPSI,'LineColor','b')
hold on
contour(0.01*TRANSP.T.PSIR,0.01*TRANSP.T.PSIZ,TPSI,'LineColor','r')
xlabel('R [m]')
ylabel('Z [m]')
title('$\psi$, blue=EFIT, red=TRANSP','Interpreter','LaTex')

% Plot gradient length of Te, should get : rise-plateau-rise-fall
dTe_drho=interpol(rho_chain2,EL.TE(ijp,:),rho_chain2,1);
figure
plot(rho_chain2,-1.*dTe_drho./EL.TE(ijp,:),'b-x')
xlabel('$\sqrt{\psi/\psi_{LCFS}}$','Interpreter','LaTex')
ylabel('$\vert \nabla(\ln(T_e)) \vert$ [m**-1]','Interpreter','LaTex')
legend('HRTS+LIDAR')


end
