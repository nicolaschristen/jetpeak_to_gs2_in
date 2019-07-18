%% Compute I(psi), required to determine Bref
%
% Input :   Rmag,Zmag -- pos. of mag. axis
%           Rflu,Zflu -- (R,Z) along chosen flux surf
%           Rpsi,Zpsi,psi -- pol. flux on fine rectangular (Rpsi,Zpsi) mesh
%           q -- safety factor on chosen flux surf.
%           psi_on_flxsurf -- value of psipol on chosen flux surf.
%           plot_verbose -- flag to produce plots to check params & fits
%
% Output:   I -- I(psi) from B = I grad(phi) + grad(phi) x grad(psi)
%
function I=compute_I(Rmag,Zmag,Rflu,Zflu,Rpsi,Zpsi,psi,q,psi_on_flxsurf,plot_verbose)

ntheta=numel(Rflu);
r=zeros(1,ntheta);
theta=zeros(1,ntheta);
dpsi_dr=zeros(1,ntheta);
dpsi_dR=zeros(1,ntheta);
dpsi_dZ=zeros(1,ntheta);

if plot_verbose
    fprintf('\nCoordinates for I(psi) integration\n')
end

for itheta=1:ntheta
    %% Compute (r,theta)
    r(itheta)=sqrt((Rflu(itheta)-Rmag)^2+ ...
        (Zflu(itheta)-Zmag)^2);
    
    % For first point, theta~0 can be slightly positive
    % or slightly negative -> use std tan function
    if itheta==1
        theta(itheta)=atan((Zflu(itheta)-Zmag)/(Rflu(itheta)-Rmag));
    % itheta=ntheta refers to the same point as itheta=1,
    % shifted by 2pi
    elseif itheta==ntheta
        theta(itheta) = theta(1) + 2*pi;
    % Otherwise, need 0<=theta<=2pi
    % -> use atan_2pi function
    else
        theta(itheta)=atan_2pi(Rflu(itheta)-Rmag, ...
            Zflu(itheta)-Zmag);
    end
    
    %% Compute dpsi_dr using chain rule
    % evaluate dpsi_dR and dpsi_dZ at closest (PSIR,PSIZ) grid-point
    deltR=Rpsi(2)-Rpsi(1);
    deltZ=Zpsi(2)-Zpsi(1);
    iR=find(abs(Rflu(itheta)-Rpsi) < (deltR/2.));
    iZ=find(abs(Zflu(itheta)-Zpsi) < (deltZ/2.));
    
    psi_iZ=psi(:,iZ); % psi at fixed Z
    psi_R_der=interpol(Rpsi,psi_iZ,Rpsi,1);
    dpsi_dR(itheta)=psi_R_der(iR);
    
    psi_iR=permute(psi,[2 1]); % psi at fixed R
    psi_iR=psi_iR(:,iR);
    psi_Z_der=interpol(Zpsi,psi_iR,Zpsi,1);
    dpsi_dZ(itheta)=psi_Z_der(iZ);
    
    dpsi_dr(itheta)=dpsi_dR(itheta)*cos(theta(itheta)) ...
        + dpsi_dZ(itheta)*sin(theta(itheta));
    
    if plot_verbose
        figure;
        plot(Rpsi,psi_iZ,'b-x')
        myxlim = xlim;
        myylim = ylim;
        hold on
        plot(Rpsi,psi_iZ(iR)+(Rpsi-Rpsi(iR))*dpsi_dR(itheta),'r-')
        xlim(myxlim)
        ylim(myylim)
        xlabel('$R$')
        ylabel('$\psi$')
        grid on
        title(['$\theta = $' num2str(theta(itheta)/pi) '$\pi$, $\Delta\psi/\psi$ = ' ...
            num2str(abs(psi_on_flxsurf-psi(iR,iZ))/psi_on_flxsurf)])
        hold on
        plot(Rpsi(iR),psi(iR,iZ),'ro')
        hold off
        
        figure;
        plot(Zpsi,psi_iR,'b-x')
        myxlim = xlim;
        myylim = ylim;
        hold on
        plot(Zpsi,psi_iR(iZ)+(Zpsi-Zpsi(iZ))*dpsi_dZ(itheta),'r-')
        xlim(myxlim)
        ylim(myylim)
        xlabel('$Z$')
        ylabel('$\psi$')
        grid on
        title(['$\theta = $' num2str(theta(itheta)/pi) '$\pi$, $\Delta\psi/\psi$ = ' ...
            num2str(abs(psi_on_flxsurf-psi(iR,iZ))/psi_on_flxsurf)])
        hold on
        plot(Zpsi(iZ),psi(iR,iZ),'ro')
        hold off
        
        fprintf('(R,Z) = (%1.3f,%1.3f)\n',Rpsi(iR),Zpsi(iZ))
    end
    
end

if plot_verbose
    fprintf('\n\n')
end

% Perform integral in theta
int_theta = r(1) / (dpsi_dr(1) * Rflu(1)) ...
        * (theta(2)-theta(1))/2.;
for itheta=2:ntheta-1
    int_theta = int_theta + ...
        r(itheta) / (dpsi_dr(itheta) * Rflu(itheta)) ...
        * (theta(itheta+1)-theta(itheta-1))/2.;
end
int_theta = int_theta + r(ntheta) / (dpsi_dr(ntheta) * Rflu(ntheta)) ...
        * (theta(ntheta)-theta(ntheta-1))/2.;

I=2.*pi*q/int_theta;


end