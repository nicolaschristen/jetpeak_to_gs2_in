function generate_gs2_infile(gs2_in,outfile_template,outfile_name)

fileID=fopen(outfile_template,'r+'); % open for reading and writing
foutID=fopen(outfile_name,'w');

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
warning(sprintf('In the original input file, parameters need to be specified in the following format, including spaces :\nname = val'))
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
            val=gs2_in.tprim1;
            if ispecies.tprim ==2
                val=gs2_in.tprim2;
            end
            ispecies.tprim=ispecies.tprim+1;
            found=1;
        case 'fprim'
            val=gs2_in.fprim1;
            if ispecies.fprim ==2
                val=gs2_in.fprim2;
            end
            ispecies.fprim=ispecies.fprim+1;
            found=1;
        case 'temp'
            val=gs2_in.temp1;
            if ispecies.temp==2
                val = gs2_in.temp2;
            end
            ispecies.temp=ispecies.temp+1;
            found=1;
        case 'dens'
            val=gs2_in.dens1;
            if ispecies.dens==2
                val = gs2_in.dens2;
            end
            ispecies.dens=ispecies.dens+1;
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