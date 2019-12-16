%% Compute the experimental radial flux given the source terms.
%
% Input :  psi -- radial profile of poloidal magnetic flux
%          dV -- radial profile of volume elements:
%                dV(i) = V(psi(i))-V(psi(i-1))
%          dx_dpsi -- dx/dpsi conversion factor to
%                     GS2 radial coordinate
%          src -- radial profile of the total source term
%
% Ouput:  flux -- experimental flux
%
function flux = flux_from_source(psi, dV, dx_dpsi, src)

nFlxSurf = numel(psi);

% Volume enclosed by each flux surface [m^3]
V = zeros(1, nFlxSurf);
V(1) = dV(1);
for iFlxSurf = 2:nFlxSurf
    V(iFlxSurf) = V(iFlxSurf-1) + dV(iFlxSurf);
end

% Compute dV/dpsi
derivOrder = 1;
dV_dpsi = interpol(psi, V, psi, derivOrder);

% Compute flux, assuming:
% src(0) = src(psimin)
% and dV/dpsi(0) = dV/dpsi(psimin)
flux = zeros(1, nFlxSurf);
flux(1) = psi(1)*dV_dpsi(1)*src(1)/2.0;
for iFlxSurf = 2:nFlxSurf
    dpsi = psi(iFlxSurf)-psi(iFlxSurf-1);
    flux(iFlxSurf) = flux(iFlxSurf-1) + ...
        dpsi * (dV_dpsi(iFlxSurf)*src(iFlxSurf) + dV_dpsi(iFlxSurf)*src(iFlxSurf))/2.0;
end
flux = 1./dV_dpsi.*flux.*dx_dpsi;

end
