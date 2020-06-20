%% Compute the experimental radial flux given the source terms.
%
% Input :  dV -- radial profile of volume elements:
%                dV(i) = V(psi(i))-V(psi(i-1))
%          dV_dpsi -- radial derivative of V
%          src -- radial profile of the total source term
%
% Ouput:  flux -- experimental flux
%
function flux = flux_from_source(dV, dV_dpsi, src)

nFlxSurf = numel(dV);

% Compute flux, assuming:
% src(0) = src(psimin)
flux = zeros(1, nFlxSurf);
flux(1) = dV(1)*src(1)/2.0;
for iFlxSurf = 2:nFlxSurf
    flux(iFlxSurf) = flux(iFlxSurf-1) + ...
        dV(iFlxSurf) * (src(iFlxSurf) + src(iFlxSurf-1))/2.0;
end
flux = 1./dV_dpsi.*flux;

end
