%% Read GS2 turbulent fluxes from a file, and give corresponding normalisations
%
% Input :   ijp --  shot index in JETPEAK DB
%           fname -- name of csv file containing computed fluxes,
%                    see example_files/fluxes.csv for format. The transport
%                    coefficients are saved in the same location as fname,
%                    in a file called 'transportCoeff.csv'
%           jData -- [kw, []] data structure obtained from JETPEAK, if
%                    read_jData has already been called.
%           trinity_norm -- [kw, 0] if true, gs2 flux dotted with gradPsi
%                           else dotted with grad(x).
%
% Ouput:    flx -- table containing fluxes read from fname
%
function flx = read_gs2Fluxes(ijp, fname, varargin)

% Read optional input arguments
options_default = struct( 'jData', [], ...
                          'trinity_norm', 0 );
opt = get_optargin(options_default, varargin);

%    ------------    %

% Read data for this shot from JETPEAK
if isempty(opt.jData)
    jData = read_jData(ijp);
else
    jData = opt.jData;
end

% Read the fluxes from file
% flux > 0 -> radially outward
flx = readtable(fname);

%    ------------    %

% Determine indices in jData corresponding to rpsi in flx
ir_jData = zeros(numel(flx.rhoc),1);
for ir_flx = 1:numel(flx.rhoc)
    r = jData.a * flx.rhoc(ir_flx);
    ir_jData(ir_flx) = find(abs(jData.rpsi-r) == min(abs(jData.rpsi-r)));
end

%    ------------    %

% Elementary charge
e=1.602e-19;

% Distinguish normalised quantities from dimensionful ones
rpsi = flx.rhoc * jData.a;
if ~opt.trinity_norm
    GammaNorm = (jData.nref(ir_jData).*jData.vthref(ir_jData))' ...
        .* jData.rhostar(ir_jData)'.^2 ...
        ./ jData.dx_dpsi(ir_jData)';
    PINorm = (jData.nref(ir_jData).*jData.mref.*jData.vthref(ir_jData).^2.*jData.a)' ...
        .* jData.rhostar(ir_jData)'.^2 ...
        ./ jData.dx_dpsi(ir_jData)';
    QNorm = (jData.nref(ir_jData)*e.*jData.tref(ir_jData) ...
            .*jData.vthref(ir_jData))' ...
        .* jData.rhostar(ir_jData)'.^2 ...
        ./ jData.dx_dpsi(ir_jData)';
else
    GammaNorm = (jData.nref(ir_jData).*jData.vthref(ir_jData))' ...
        .* jData.rhostar(ir_jData)'.^2 ...
        .* jData.gradPsiAvg(ir_jData)';
    PINorm = (jData.nref(ir_jData).*jData.mref.*jData.vthref(ir_jData).^2.*jData.a)' ...
        .* jData.rhostar(ir_jData)'.^2 ...
        .* jData.gradPsiAvg(ir_jData)';
    QNorm = (jData.nref(ir_jData)*e.*jData.tref(ir_jData) ...
            .*jData.vthref(ir_jData))' ...
        .* jData.rhostar(ir_jData)'.^2 ...
        .* jData.gradPsiAvg(ir_jData)';
end

%    ------------    %

% Append information for normalisation to table
to_add_to_table = table(rpsi,GammaNorm,PINorm,QNorm,ir_jData);
flx = [flx to_add_to_table];

end
