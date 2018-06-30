%% Compute I(psi), required to determine Bref
% Input :   Rmag,Zmag -- pos. of mag. axis
%           Rflu,Zflu -- (R,Z) along chosen flux surf
%           Rpsi,Zpsi,psi -- pol. flux on fine rectangular (Rpsi,Zpsi) mesh
%           q -- safety factor on chosen flux surf.
function I=compute_I(Rmag,Zmag,Rflu,Zflu,Rpsi,Zpsi,psi,q)

ntheta=numel(Rflu);
r=zeros(1,ntheta);
theta=zeros(1,ntheta);
dpsi_dr=zeros(1,ntheta);
dpsi_dR=zeros(1,ntheta);
dpsi_dZ=zeros(1,ntheta);
int_theta=0.;

for itheta=1:ntheta
    %% Compute jacobian Jr from (r,theta,zeta) to (R,Z,zeta)
    r(itheta)=sqrt((Rflu(itheta)-Rmag)^2+ ...
        (Zflu(itheta)-Zmag)^2);
    
    theta(itheta)=atan_2pi(Rflu(itheta)-Rmag, ...
        Zflu(itheta)-Zmag);
    % Correct last element of theta
    % because atan_2pi returns theta(1)=theta(ntheta)
    if itheta==ntheta 
        theta(itheta)=2.*pi+theta(1);
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
    
end

% Perform integral in theta
int_theta = int_theta + r(1) / (dpsi_dr(1) * Rflu(1)) ...
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