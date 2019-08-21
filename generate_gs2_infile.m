%% Modify GS2 template input file
%
% Input :   gs2_in -- structure with dimensionless GS2 parameters
%           infile_template -- GS2 input file template
%           new_infile_name -- name of input file to be generated
%           add_carbon -- =0 neglects C and adapts n_ion so that QN holds
%                         =1 adds Carbon as an impurity
%
% Output:   -
%
function generate_gs2_infile(gs2_in,infile_template,new_infile_name,add_carbon)

fileID=fopen(infile_template,'r+'); % open for reading and writing
foutID=fopen(new_infile_name,'w');

iline=0;
line=fgetl(fileID);
while ischar(line)
    content{iline+1}=line;
    iline=iline+1;
    line=fgetl(fileID);
end
nline=iline;

%celldisp(content)

found=0;
ispecies.tprim=1;
ispecies.fprim=1;
ispecies.dens=1;
ispecies.temp=1;
ispecies.vnewk=1;
fprintf(['In the template input file, parameters need to be specified \n', ...
    'in the following format (including spaces) :\n', ...
    'name = val\n\n'])
if add_carbon
    fprintf(['Carbon is being included. Please make sure that \n'...
        'the template input file has 3 sets of species namelists'])
end
for iline=1:nline
    try
        [name, val]=strread(content{iline},'%s = %f');
    switch name{1}
        case 'Rmaj'
            val = gs2_in.Rmaj;
            found=1;
        case 'R_geo'
            val = gs2_in.R_geo;
            found=1;
        case 'rhoc'
            val = gs2_in.rhoc;
            found=1;
        case 'shift'
            val = gs2_in.shift;
            found=1;
        case 'qinp'
            val = gs2_in.qinp;
            found=1;
        case 'shat'
            val = gs2_in.shat;
            found=1;
        case 's_hat_input'
            val = gs2_in.shat;
            found=1;
        case 'akappa'
            val = gs2_in.akappa;
            found=1;
        case 'akappri'
            val = gs2_in.akappri;
            found=1;
        case 'tri'
            val = gs2_in.tri;
            found=1;
        case 'tripri'
            val = gs2_in.tripri;
            found=1;
        case 'beta'
            val = gs2_in.beta;
            found=1;
        case 'zeff'
            val = gs2_in.zeff;
            found=1;
        case 'beta_prime_input'
            val = gs2_in.beta_prime_input;
            found=1;
        case 'mach'
            val = gs2_in.mach;
            found=1;
        case 'g_exb'
            val = gs2_in.g_exb;
            found=1;
        case 'tprim'
            switch ispecies.tprim
                case 1
                    val=gs2_in.tprim1;
                case 2
                    val=gs2_in.tprim2;
                case 3
                    val=gs2_in.tprim3;
            end
            ispecies.tprim=ispecies.tprim+1;
            found=1;
        case 'fprim'
            switch ispecies.fprim
                case 1
                    val=gs2_in.fprim1;
                case 2
                    val=gs2_in.fprim2;
                case 3
                    val=gs2_in.fprim3;
            end
            ispecies.fprim=ispecies.fprim+1;
            found=1;
        case 'temp'
            switch ispecies.temp
                case 1
                    val=gs2_in.temp1;
                case 2
                    val=gs2_in.temp2;
                case 3
                    val=gs2_in.temp3;
            end
            ispecies.temp=ispecies.temp+1;
            found=1;
        case 'dens'
            switch ispecies.dens
                case 1
                    val=gs2_in.dens1;
                case 2
                    val=gs2_in.dens2;
                case 3
                    val=gs2_in.dens3;
            end
            ispecies.dens=ispecies.dens+1;
            found=1;
        case 'vnewk'
            switch ispecies.vnewk
                case 1
                    val=gs2_in.vnewk1;
                case 2
                    val=gs2_in.vnewk2;
                case 3
                    val=gs2_in.vnewk3;
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

end