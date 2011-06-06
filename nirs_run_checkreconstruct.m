function out = nirs_run_checkreconstruct(job)
%
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire
% Clement Bonnery
Kas =3;
switch Kas
    case 3
        out_Hbs = job.outreconstruct_Hb;
        load(job.NIRSmat{:});
        
        for iHbs =1:size(out_Hbs,1)
            f = out_Hbs{iHbs,1};
            V = spm_vol(f);
            Y = spm_read_vols(V);
            
            if iHbs==1
                [dirr,dummy,dummy2] = fileparts(f);
                Y_t_O = zeros([size(Y) size(out_Hbs,1)/2]);
                Y_t_R = zeros([size(Y) size(out_Hbs,1)/2]);
            end
            
            [dummy,name,dummy2] = fileparts(f);
            sep =strfind(name,'_');
            %             timee =str2num(name(sep(2)+2:sep(3)-1));
            timee =str2num(name(sep(2)+2:end));
            itO =1;
            itR =1;
            if ~isempty(strfind(name,'HbO'))
                Y_t_O(:,:,:,itO) = Y;
                itO = itO+1;
            else
                Y_t_R(:,:,:,itR) = Y;
                itR = itR+1;
            end 
        end
        
       %calcul de la variance sur chaque voxel pour chaque condition
       sigmaO_baseline = zeros(size(Y));
       sigmaR_baseline = zeros(size(Y));
       
       sigmaO_conditi1 = zeros(size(Y));
       sigmaR_conditi1 = zeros(size(Y));
       
       XO_baseline = zeros(size(Y));
       XR_baseline = zeros(size(Y));
       
       XO_conditi1 = zeros(size(Y));
       XR_conditi1 = zeros(size(Y));
       
       for i=1:size(Y,1)
           for j=1:size(Y,2)
               for k=1:size(Y,3)
                   %
                   sigmaO_baseline(i,j,k) = std(Y_t_O(i,j,k,1:66));
                   XO_baseline(i,j,k) = mean(Y_t_O(i,j,k,1:66));
                   
                   sigmaR_baseline(i,j,k) = std(Y_t_R(i,j,k,1:66));
                   XR_baseline(i,j,k) = mean(Y_t_R(i,j,k,1:66));
                   
                   %
                   sigmaO_conditi1(i,j,k) = std(Y_t_O(i,j,k,67:132));
                   XO_conditi1(i,j,k) = mean(Y_t_O(i,j,k,67:132));
                   
                   sigmaR_conditi1(i,j,k) = std(Y_t_R(i,j,k,67:132));
                   XR_conditi1(i,j,k) = mean(Y_t_R(i,j,k,67:132));
               end
           end
       end
       
       for i=1:size(sigmaO_baseline,1)
       sigmaO(i,:,:) = sqrt(1/2*(squeeze(sigmaO_baseline(i,:,:)).^2+squeeze(sigmaO_conditi1(i,:,:)).^2));
       sigmaR(i,:,:) = sqrt(1/2*(squeeze(sigmaR_baseline(i,:,:)).^2+squeeze(sigmaR_conditi1(i,:,:)).^2));
       end
       
       Y_st_O = (XO_conditi1 - XO_baseline)./sigmaO;
       Y_st_R = (XR_conditi1 - XR_baseline)./sigmaR;
       
       Y_st_O(isnan(Y_st_O))=0;
       Y_st_R(isnan(Y_st_R))=0;
       
        %%% nifti
        V_O = struct('fname',fullfile(dirr,'stat_t_HbO.nii'),...
            'dim',  V.dim,...
            'dt',   V.dt,...
            'pinfo',V.pinfo,...
            'mat',  V.mat);
        
        V_O = spm_create_vol(V_O);
        spm_write_vol(V_O, Y_st_O);
        
        V_R = struct('fname',fullfile(dirr,'stat_t_HbR.nii'),...
            'dim',  V.dim,...
            'dt',   V.dt,...
            'pinfo',V.pinfo,...
            'mat',  V.mat);

        V_R = spm_create_vol(V_R);
        spm_write_vol(V_R, Y_st_R);
        
    case 1%%%%% variations temporelles sur un voxel
        out_Hbs = job.outreconstruct_Hb;
        
        load(job.NIRSmat{:});
        fnirs = load(NIRS.Dt.fir.pp.p{:},'-mat');
        tf = length(fnirs.d(:,1));
        t = (1:tf)*0.04;
        
        dec_DHbO = zeros(1,length(fnirs.d));
        dec_DHbR = zeros(1,length(fnirs.d));
        
        for iHbs =1:size(out_Hbs,1)
            f = out_Hbs{iHbs,1};
            V = spm_vol(f);
            Y = spm_read_vols(V);
            
            [dummy,name,dummy2] = fileparts(f);
            sep =strfind(name,'_');
            %             timee =str2num(name(sep(2)+2:sep(3)-1));
            timee =str2num(name(sep(2)+2:end));
            if ~isempty(strfind(name,'HbO'))
                
                dec_DHbO(1,timee:timee+25) = Y(13,7,8);
            else
                dec_DHbR(1,timee:timee+25) = Y(13,7,8);
            end
        end
        figure;
        subplot(2,1,1)
        plot(t,dec_DHbO,'r');
        subplot(2,1,2)
        plot(t,dec_DHbR);
        
    case 2
        out_muas = job.outreconstruct_Hb;
        
        load(job.NIRSmat{:});
        fnirs = load(NIRS.Dt.fir.pp.p{:},'-mat');
        
        dec_Dmua = zeros(1,length(fnirs.d));
        
        for imuas =1:size(out_muas,1)
            f = out_muas{imuas,1};
            V = spm_vol(f);
            Y = spm_read_vols(V);
            
            [dummy,name,dummy2] = fileparts(f);
            sep =strfind(name,'_');
            timee =str2num(name(sep(2)+2:sep(3)-1));
            
            dec_Dmua(1,timee:timee+25) = Y(11,8,7);
        end
        figure;
        plot(dec_Dmua(1,5000:7000));
end
out =1;
end