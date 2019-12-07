%% Fill GS2 template input file with gs2 variables extracted from JETPEAK/TRANSP
%
% Input :   ijp -- JETPEAK index of shot
%           rpsi_nrm -- approx. rpsi/a value(s) at which the file should be generated
%                       (will choose nearest grid points, no interpolation)
%           template_fName -- name of GS2 input file template
%           new_fName -- name(s) of input file to be generated (same size as rpsi_nrm)
%                        if multiple names, pass them as a cell-array.
%
% Output:   -
%
function create_gs2VarsFile(ijp,rpsi_nrm,template_fName,new_fName)

if ~iscell(new_fName)
    new_fName = {new_fName};
end

if numel(rpsi_nrm) ~= numel(new_fName)
    disp(['Mismatch between the number of radial locations,\n', ...
        'and the number of new file names. Exiting.'])
    exit
end

nFiles = numel(rpsi_nrm);

% Extract profiles from JETPEAK/TRANSP
jData = read_jData(ijp);

% Compute profiles of gs2 variables from this data
gs2Vars = gs2Vars_from_jData(jData);

for iFile = 1:nFiles

    fileID=fopen(template_fName,'r+'); % open for reading and writing
    foutID=fopen(new_fName{iFile},'w');

    % Find flux surface nearest to the user's choice
    [~, iFlxSurf] = min(abs(gs2Vars.rhoc - rpsi_nrm(iFile)));
    
    iline=0;
    line=fgetl(fileID);
    while ischar(line)
        content{iline+1}=line;
        iline=iline+1;
        line=fgetl(fileID);
    end
    nline=iline;
    
    found=0;
    ispecies.tprim=1;
    ispecies.fprim=1;
    ispecies.dens=1;
    ispecies.temp=1;
    ispecies.vnewk=1;
    fprintf(['In the template input file, parameters need to be specified \n', ...
        'in the following format (including spaces) :\n', ...
        'name = val\n\n'])
    if gs2Vars.add_carbon(iFlxSurf)
        fprintf(['Carbon is being included. Please make sure that \n'...
            'the template input file has 3 sets of species namelists'])
    end
    for iline=1:nline
        try
            [name, val]=strread(content{iline},'%s = %f');
        switch name{1}
            case 'Rmaj'
                val = gs2Vars.Rmaj(iFlxSurf);
                found=1;
            case 'R_geo'
                val = gs2Vars.R_geo(iFlxSurf);
                found=1;
            case 'rhoc'
                val = gs2Vars.rhoc(iFlxSurf);
                found=1;
            case 'shift'
                val = gs2Vars.shift(iFlxSurf);
                found=1;
            case 'qinp'
                val = gs2Vars.qinp(iFlxSurf);
                found=1;
            case 'shat'
                val = gs2Vars.shat(iFlxSurf);
                found=1;
            case 's_hat_input'
                val = gs2Vars.shat(iFlxSurf);
                found=1;
            case 'akappa'
                val = gs2Vars.akappa(iFlxSurf);
                found=1;
            case 'akappri'
                val = gs2Vars.akappri(iFlxSurf);
                found=1;
            case 'tri'
                val = gs2Vars.tri(iFlxSurf);
                found=1;
            case 'tripri'
                val = gs2Vars.tripri(iFlxSurf);
                found=1;
            case 'beta'
                val = gs2Vars.beta(iFlxSurf);
                found=1;
            case 'zeff'
                val = gs2Vars.zeff(iFlxSurf);
                found=1;
            case 'beta_prime_input'
                val = gs2Vars.beta_prime_input(iFlxSurf);
                found=1;
            case 'mach'
                val = gs2Vars.mach(iFlxSurf);
                found=1;
            case 'g_exb'
                val = gs2Vars.g_exb(iFlxSurf);
                found=1;
            case 'tprim'
                switch ispecies.tprim
                    case 1
                        val=gs2Vars.tprim1(iFlxSurf);
                    case 2
                        val=gs2Vars.tprim2(iFlxSurf);
                    case 3
                        val=gs2Vars.tprim3(iFlxSurf);
                end
                ispecies.tprim=ispecies.tprim+1;
                found=1;
            case 'fprim'
                switch ispecies.fprim
                    case 1
                        val=gs2Vars.fprim1(iFlxSurf);
                    case 2
                        val=gs2Vars.fprim2(iFlxSurf);
                    case 3
                        val=gs2Vars.fprim3(iFlxSurf);
                end
                ispecies.fprim=ispecies.fprim+1;
                found=1;
            case 'temp'
                switch ispecies.temp
                    case 1
                        val=gs2Vars.temp1(iFlxSurf);
                    case 2
                        val=gs2Vars.temp2(iFlxSurf);
                    case 3
                        val=gs2Vars.temp3(iFlxSurf);
                end
                ispecies.temp=ispecies.temp+1;
                found=1;
            case 'dens'
                switch ispecies.dens
                    case 1
                        val=gs2Vars.dens1(iFlxSurf);
                    case 2
                        val=gs2Vars.dens2(iFlxSurf);
                    case 3
                        val=gs2Vars.dens3(iFlxSurf);
                end
                ispecies.dens=ispecies.dens+1;
                found=1;
            case 'vnewk'
                switch ispecies.vnewk
                    case 1
                        val=gs2Vars.vnewk1(iFlxSurf);
                    case 2
                        val=gs2Vars.vnewk2(iFlxSurf);
                    case 3
                        val=gs2Vars.vnewk3(iFlxSurf);
                end
                ispecies.vnewk=ispecies.vnewk+1;
                found=1;
        end
        catch ME
            if strcmp(ME.identifier,'MATLAB:dataread:TroubleReading')
                % Line with different formatting than expected -> don't edit
                found=0;
            elseif strcmp(ME.identifier,'MATLAB:badsubscript')
                % Empty line -> don't edit
                found=0;
            end
        end
        if found==1
            content{iline}=[' ' name{1} ' = ' num2str(val)];
        end
        found=0;
        fprintf(foutID,'%s\r\n',content{iline});
    end
    
    fclose(foutID);
    fclose(fileID);

end % loop over files

end
