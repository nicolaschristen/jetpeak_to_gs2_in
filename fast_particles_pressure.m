load databases/TRANSP_2017_3.mat
load databases/JETPEAK_2019_10.mat

ijp = 950;
idxAscot = find(ASCOT.SAMPLE==ijp); % might find multiple instances
idxTRANSP = find(TRANSP.JPI==ijp);

jData = read_jData(ijp);

e = 1.6e-19; % elementary charge
rho_chain2=linspace(0.,1.,21); % sqrt(psin)

% ASCOT
p_fast_ASC = ASCOT.PRFAST(idxAscot,:);
p_fast_ASC = squeeze(mean(p_fast_ASC,1)); % average over different ASCOT instances
p_fast_ASC = interpol(rho_chain2, p_fast_ASC, jData.rpsi); % interpolate

% TRANSP
p_fast_TRANSP = 2/3 * squeeze(TRANSP.F.UFIPP(idxTRANSP,:) ...
    + TRANSP.F.UFIPA(idxTRANSP,:)); % sum of parallel and perp energy dens
TRANSP_has_fast = any(p_fast_TRANSP ~= 0);

figure
plot(jData.rpsi,jData.ni.*jData.ti*e)
hold on
plot(jData.rpsi,p_fast_ASC)
if TRANSP_has_fast
    hold on
    plot(jData.rpsi,p_fast_TRANSP)
end
grid on
xlabel('$r_\psi$ [m]')
ylabel('$p$ [Pa]')
if TRANSP_has_fast
    legend('Thermal', 'Fast (ASCOT)', 'Fast (TRANSP)')
else
    legend('Thermal', 'Fast (ASCOT)')
end

