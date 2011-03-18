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
cs = NIRS.Cs.mcs{ics};

if cs.alg==1
    Oe='.inp';
elseif cs.alg==2
    Oe='.2pt';
end
segR = cs.segR;
V_segR = spm_vol(segR);

%%%% lecture du bin pour les Boas parts
b8i = cs.b8i;
% V_b8i = spm_vol(b8i);
%%%%%% trouver LA TAILLE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%%% b8i a des valeurs NULLES !!!!!!!!!!!!!!!!!!!!!!!!
%%%%%% VALEURS DES MD ET MS !!!!!!!!!!!!!!!!!!!!!!!!!!!!
fid=fopen(b8i,'rb');
b8i = fread(fid,11*12*13,'double');
fclose(fid);
% b8i = reshape(b8i,V_b8i);
%%%%

wl = NIRS.Cf.dev.wl; %= [830 690];
Cid = NIRS.Cf.H.C.id;
Cwl = NIRS.Cf.H.C.wl;
C =[];
ND =1;
isd =0; % count of lignes in sensitivity matrix

for fi = 1:size(f,1)
    [~,fn,~] = fileparts(f{fi,:});
    [~,Pwl] = find(wl==str2double(fn(8:10)));
    
    if strcmp(fn(1:1),'D')% stores all ''D_No.._...nm.2pt'' files in a list
        fD{ND,:} = f{fi,:};
        ND = ND+1;
        
    elseif strcmp(fn(1:1),'S')
        Sn = str2double(fn(5:6));
        
        fid=fopen(f{fi,:},'rb');
        % in Boas code : float 64 : but same result...
        ms = fread(fid,V_segR.dim(1)*V_segR.dim(2)*V_segR.dim(3),'double');
        fclose(fid);
        % Boas : scaling
%         ms = reshape(ms,V_segR.dim);
        ms = ms / cs.par.nphotons;
        lst = find(ms<0); ms_neg = sum(ms(lst));
        lst = find(ms>0); ms_pos = sum(ms(lst).*b8i(lst));
        ms(lst) = ms(lst) * (1+ms_neg) / ms_pos;
        % Boas : end scaling
       
        D_Sn =  Cid(3,Cid(2,:)== Sn);% Detectors seeing source Sn
        for i = 1:size(D_Sn,2)% overview of detectors seen by source thanks to Cid
             for j = 1:size(fD,1)% only if detector has been selected by user...
                [~,fD_n,~] = fileparts(fD{j,:});
                
                if D_Sn(1,i) < 10, Ds_Sn = ['0' int2str(D_Sn(1,i))]; else Ds_Sn = int2str(D_Sn(1,i));end
                Dfn = ['D_No' Ds_Sn '_' int2str(wl(Pwl)) 'nm'];
                b = strcmp(fD_n,Dfn);
                if b % one pair found : processing goes on !
                    c = find(Cid(2,:)== Sn & Cid(3,:)==D_Sn(i) & Cwl==Pwl);
                    
                    fid=fopen(fullfile(cs_dir,[Dfn Oe]),'rb');
                    md = fread(fid,V_segR.dim(1)*V_segR.dim(2)*V_segR.dim(3),'double');
                    fclose(fid);
                    % Boas : scaling
                    md = md / cs.par.nphotons;
                    lst = find(ms<0); md_neg = sum(md(lst));
                    lst = find(ms>0); md_pos = sum(md(lst).*b8i(lst));
                    md(lst) = md(lst) * (1+md_neg) / md_pos;
                    
                    phi0 = md( floor(pos(idxD,1)), floor(pos(idxD,2)),floor(pos(idxD,3)) );
                end
                % Boas : end scaling
                
                sens_sd = ms.*md;
                
                for idx=1:nMeas
                    A(:,:,:,idx) = p2pt(:,:,:,ml(idx,1)) .* p2pt(:,:,:,ml(idx,2)) ...
                        / ((phio(ml(idx,1),ml(idx,2)) + phio(ml(idx,2),ml(idx,1))) / 2);
                end
                
                sens(isd+1,:) = sens_sd';%double();
                isd = isd+1;
                %%%%% on precise les paires qu'on a
                C = [C c];
            end
        end
    end
end
save(fullfile(cs_dir,'sens.mat'),'sens');

cs.C = C;
save(job.NIRSmat{1,1},'NIRS');

%%% au passage on genere un nii, a mettre en OPTION...
%%% RESHAPE_SENSITIVITY_MATRIX
sens_reshaped = zeros(V_segR.dim);
for i=1:size(sens,1)
    sens_reshaped = sens_reshaped + reshape(sens(i,:),V_segR.dim);
end

V = struct('fname',fullfile(cs_dir,['sens' '.nii']),...
    'dim',  V_segR.dim,...
    'dt',   V_segR.dt,...
    'pinfo',V_segR.pinfo,...
    'mat',  V_segR.mat);

V = spm_create_vol(V);
V = spm_write_vol(V, sens_reshaped);
%%%

out = 1;
end