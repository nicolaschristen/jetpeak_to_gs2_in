%% Compute the flux-surface averaged norm of the gradient
% of the poloidal magnetic flux.
%
% Input : Rpsi -- major radius for rectangular grid
%         Zpsi -- vertical coord. for rectangular grid
%         psi -- poloidal magnetic flux on rectangular grid
%         Rflu -- major radius on (psi,theta) grid
%         Zflu -- vertical coord. on (psi,theta) grid
%         Rmag -- major radius of the magnetic axis
%         Zmag -- vertical position of the magnetic axis
%
% Output: gradPsiAvg -- flux-surface averaged norm of the gradient
%                       of the poloidal magnetic flux, as a function
%                       of radial coordinate.
%
function gradPsiAvg = compute_gradPsiAvg(Rpsi, Zpsi, psi, Rflu, Zflu, Rmag, Zmag)

nR = numel(Rpsi);
dR = Rpsi(2)-Rpsi(1);
nZ = numel(Zpsi);
dZ = Zpsi(2)-Zpsi(1);
nflxsurf = size(Rflu,1);
ntheta = size(Rflu,2);

% Compute dpsi/dR on rectangular mesh
% with centered differences, except for edges of mesh
dpsi_dR = zeros(nR,nZ);
for iZ = 1:nZ
    dpsi_dR(1,iZ) = (psi(2,iZ)-psi(1,iZ))/dR;
    for iR = 2:nR-1
        dpsi_dR(iR,iZ) = (psi(iR+1,iZ)-psi(iR-1,iZ))/(2*dR);
    end
    dpsi_dR(nR,iZ) = (psi(nR,iZ)-psi(nR-1,iZ))/dR;
end

% Compute dpsi/dZ on rectangular mesh
% with centered differences, except for edges of mesh
dpsi_dZ = zeros(nR,nZ);
for iR = 1:nR
    dpsi_dZ(iR,1) = (psi(iR,2)-psi(iR,1))/dZ;
    for iZ = 2:nZ-1
        dpsi_dZ(iR,iZ) = (psi(iR,iZ+1)-psi(iR,iZ-1))/(2*dZ);
    end
    dpsi_dZ(iR,nZ) = (psi(iR,nZ)-psi(iR,nZ-1))/dZ;
end

% Compute theta
theta = zeros(nflxsurf,ntheta);
for iflx = 1:nflxsurf
    for itheta=1:ntheta
        % For first point, theta~0 can be slightly positive
        % or slightly negative -> use std tan function
        if itheta==1
            theta(iflx,itheta)=atan((Zflu(iflx,itheta)-Zmag)/(Rflu(iflx,itheta)-Rmag));
        % itheta=ntheta refers to the same point as itheta=1,
        % shifted by 2pi
        elseif itheta==ntheta
            theta(iflx,itheta) = theta(iflx,1) + 2*pi;
        % Otherwise, need 0<=theta<=2pi
        % -> use atan_2pi function
        else
            theta(iflx,itheta)=atan_2pi(Rflu(iflx,itheta)-Rmag, ...
                Zflu(iflx,itheta)-Zmag);
        end
    end
end

% Compute |gradPsi| and average over flux-surfaces
gradPsi = @(iR,iZ) sqrt(dpsi_dR(iR,iZ)^2 + dpsi_dZ(iR,iZ)^2);
gradPsiAvg = zeros(1,nflxsurf);
get_near_iR = @(iflx,itheta) find(abs(Rflu(iflx,itheta)-Rpsi)<dR/2);
get_near_iZ = @(iflx,itheta) find(abs(Zflu(iflx,itheta)-Zpsi)<dZ/2);
for iflx = 1:nflxsurf
    iRnear = get_near_iR(iflx,1);
    iZnear = get_near_iZ(iflx,1);
    intgrl = gradPsi(iRnear,iZnear) ...
        * (theta(iflx,2)-theta(iflx,1))/2.;
    for itheta = 2:ntheta-1
        iRnear = get_near_iR(iflx,itheta);
        iZnear = get_near_iZ(iflx,itheta);
        intgrl = intgrl + gradPsi(iRnear,iZnear) ...
            * (theta(iflx,itheta+1)-theta(iflx,itheta-1))/2.;
    end
    iRnear = get_near_iR(iflx,ntheta);
    iZnear = get_near_iZ(iflx,ntheta);
    intgrl = intgrl + gradPsi(iRnear,iZnear) ...
        * (theta(iflx,ntheta)-theta(iflx,ntheta-1))/2.;
    gradPsiAvg(iflx) = intgrl/(2.*pi);
end

end
