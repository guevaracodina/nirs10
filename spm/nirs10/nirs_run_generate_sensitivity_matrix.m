function out = nirs_run_generate_sensitivity_matrix(job)
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________

%Usage: only works with .mc2 or .2pt files for now

%Load NIRS.mat information
load(job.NIRSmat{1,1});

f = job.outMCfiles;
cs_dir =  fileparts(f{1,:});
cs_ldir = cs_dir(max(strfind(cs_dir,'\'))+9:end);

ics =1;
while ~strcmp(cs_ldir,NIRS.Cs.n{ics}) 
    ics =ics+1; 
end

if NIRS.Cs.mcs{ics}.alg==1
    Oe='.'; 
elseif NIRS.Cs.mcs{ics}.alg==2
    Oe='.2pt'; 
end
segR = NIRS.Cs.mcs{ics}.segR;
V_segR = spm_vol(segR);

wl = NIRS.Cf.dev.wl; %= [830 690];
Cid = NIRS.Cf.H.C.id;
Cwl = NIRS.Cf.H.C.wl;

count =1;
isd =0; % number of lignes in sensitivity matrix

for fi = 1:size(f,1)
    [~,fn,~] = fileparts(f{fi,:});
    [~,Pwl] = find(wl==str2double(fn(8:10)));
    
    if strcmp(fn(1:1),'D')% tous les 2pt des detecteurs sont ranges parce que D precede S dans l'alphabet
        fD{count,:} = f{fi,:};
        count = count+1;
        
    elseif strcmp(fn(1:1),'S')
        Sn_fi = str2double(fn(5:6));

        fid=fopen(f{fi,:},'rb');
        ms = fread(fid,V_segR.dim(1)*V_segR.dim(2)*V_segR.dim(3),'double');%'single');
        fclose(fid);
        
        [~,D_Sn] = find(Cid(2,:)== Sn_fi);% looking for detectors seeing current source
        for i = 1:size(D_Sn,2)% overview of detectors seen by source thanks to Cid
            for j = 1:size(fD,1)% only if detector has been selected by user...
                [~,fD_n,~] = fileparts(fD{j,:});
                
                if D_Sn(1,i) < 10, Ds_Sn = ['0' int2str(D_Sn(1,i))]; else Ds_Sn = int2str(D_Sn(1,i));end
                Dfn = ['D_No' Ds_Sn '_' int2str(wl(Pwl)) 'nm'];
                b = strcmp(fD_n,Dfn);
                if b % one pair found : processing goes on !
                    [~,wl_Sn] = find(Cwl==Pwl);
                    [~,c] = find(wl_Sn==D_Sn(1,i));
%                     Cn_it = wl_Sn(c);
                    
                    fid=fopen(fullfile(cs_dir,[Dfn Oe]),'rb');
                    md = fread(fid,V_segR.dim(1)*V_segR.dim(2)*V_segR.dim(3),'double');%'single');
                    fclose(fid);
                    
                    %sens_sd = real(ms)/max(real(ms)).*real(md)/max(real(md)); %treated as long vectors
                    sens_sd = ms.*md;
                    sens_sd = log(sens_sd);
                    
                    sens(isd+1,:) = sens_sd';%double();
                    isd = isd+1;
                end
            end
        end
    end
end
save(fullfile(cs_dir,'sens.mat'),'sens');

out = 1;
end

%         %meme caracteristiques que l'image binaire envoyee
%         cs = NIRS.Cs.mcs{ics,1};
%         [~,n_segR,~] = fileparts(cs.segR);
%         V_segR = spm_vol(cs.segR);
%
%         sens_sd = reshape(sens_sd,V_segR.dim);
%
%         V_sens.fname = fullfile(cs_dir,['sens' n_segR(1:11) '.nii']);
%         V_sens.dim   = V_segR.dim; %same as dimi
%         V_sens.dt    = V_segR.dt; %float32 not V.dt which is int16
%         V_sens.mat   = V_segR.mat;
%         V_sens.pinfo = V_segR.pinfo;
%         spm_write_vol(V_sens,sens_sd);
%
%         jobR.image_in{:} = V_sens.fname;
%         jobR.out_dt = 'same';
%         jobR.out_autonaming = 0;
%         jobR.out_prefix = 'bla';
%         jobR.out_dim = V_seg.dim;
%         jobR.out_dir = cs_dir;
%         outR = nirs_resize(jobR);
%
%         %                     V.fname = [dir1 '\srcNo' file_str_s '_' num2str(wl(iwl)) 'nm.nii'];
%         %                     ms_temp = reshape(log(ms),dimf);
%         %                     ms_temp = interp3(xf1,yf1,zf1,ms_temp,xi1,yi1,zi1,'nearest');
%         %                     spm_write_vol(V,ms_temp);
%         %                     V.fname = [dir1 '\detNo' file_str_s '_' num2str(wl(iwl)) 'nm.nii'];
%         %                     md_temp = reshape(log(md),dimf);
%         %                     md_temp = interp3(xf1,yf1,zf1,md_temp,xi1,yi1,zi1,'nearest');
%         %                     spm_write_vol(V,md_temp);
%         %
%         %Save in .nii format - use structure of segmented file
%         V.fname = filen_sd; %[dir1 filen_sd];
%         V.dim  = V.dim; %same as dimi
%         V.dt   = [16 0]; %float32 not V.dt which is int16
%         V.mat   = V.mat;
%         V.pinfo = V.pinfo;
%         spm_write_vol(V,outR);
