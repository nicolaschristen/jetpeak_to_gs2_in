%% Compute Coulomb logarithm for electron-ion collisions
%
% Input :   Zi -- ion atomic number
%           mi -- ion mass [kg]
%           ne -- radial profile of electron density [m^{-3}]
%           ni -- radial profile of ion density [m^{-3}]
%           te -- radial profile of electron temperature [eV]
%           ti -- radial profile of ion temperature [eV]
%
% Output:   loglambda -- radial profile of the Coulomb logarithm
%
function loglambda = coulomb_log(Zi,mi,ne,ni,te,ti)

% Physical constants
mp = 1.673e-27; % proton mass
me = 9.109e-31; % electron mass

% Unit conversion
cm_per_m = 10^2;

Nrho = numel(ne); % Number of radial points in a profile
loglambda = 17*ones(1,Nrho);

for irho=1:Nrho
    % Convert dens to [cm^{-3}] (see NRL formulary)
    my_ne = ne(irho) * cm_per_m^(-3);
    my_ni = ni(irho) * cm_per_m^(-3);
    my_te = te(irho);
    my_ti = ti(irho);
    
    if my_ti*me/mi < my_te && my_te < 10*Zi^2
        loglambda(irho) = 23 - log(my_ne^0.5*Zi*my_te^(-1.5));
    elseif my_ti*me/mi < 10*Zi^2 && 10*Zi^2 < my_te
        loglambda(irho) = 24 - log(my_ne^0.5*my_te^(-1));
    elseif my_te < my_ti*Zi*me/mi
        loglambda(irho) = 30 - log(my_ni^0.5*my_ti^(-1.5)*Zi^2*mp/mi);
    else
        fprintf('WARNING in coulomb_log, index %u: this case is not covered.\n',irho)
        fprintf('\tSetting loglambda(%u) = 17\n\n',irho)
    end
    
end

end