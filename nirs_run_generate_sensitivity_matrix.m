function out = nirs_run_generate_sensitivity_matrix(job)
%
%
% Clément Bonnéry 03/2011

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
    Oe='.mc2';
elseif cs.alg==2
    Oe='.2pt';
end
segR = cs.segR;
V_segR = spm_vol(segR);

% lecture du fichier binaire pour avoir les mu_a
b8i = cs.b8i;
fid=fopen(b8i,'rb');
Yb8i = fread(fid,prod(V_segR.dim),'uint8');
fclose(fid);
% : mua mus' ... on prend les mu_a
Yb8i_l1 = zeros(size(Yb8i));
Yb8i_l1(Yb8i==1)= cs.par.gmPpties_l1(1,1);
Yb8i_l1(Yb8i==2)= cs.par.wmPpties_l1(1,1);
Yb8i_l1(Yb8i==3)= cs.par.csfPpties_l1(1,1);
Yb8i_l1(Yb8i==4)= cs.par.skullPpties_l1(1,1);
Yb8i_l1(Yb8i==5)= cs.par.scalpPpties_l1(1,1);

Yb8i_l2 = zeros(size(Yb8i));
Yb8i_l2(Yb8i==1)= cs.par.gmPpties_l2(1,1);
Yb8i_l2(Yb8i==2)= cs.par.wmPpties_l2(1,1);
Yb8i_l2(Yb8i==3)= cs.par.csfPpties_l2(1,1);
Yb8i_l2(Yb8i==4)= cs.par.skullPpties_l2(1,1);
Yb8i_l2(Yb8i==5)= cs.par.scalpPpties_l2(1,1);

if size(Yb8i_l1,1)==1%cs.alg==2
    Yb8i_l1 = Yb8i_l1';
    Yb8i_l2 = Yb8i_l2';
end

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
        
        switch cs.alg
            case 1 % MCX
                ms=loadmc2(f{fi,:},V_segR.dim);
                % pour MCX, on obtient P(r,t) qu'il faut nomaliser avec la
                % formule (1) de l'article de Fang (si on lit
                %http://mcx.sourceforge.net/cgi-bin/index.cgi?MMC/Doc/FAQ#How_do_I_interpret_MMC_s_output_data
                % , on comprends )
                % Bonnery from Boas et al.
                %                 sum_Jout = sum(ms(ms<0)/cs.par.nphotons);
                %                 Svx = 5*5;
                %                 Vvx = 5*5*5;
                %                 msP = ms;
                %                 msP(msP<0)=0;
                %                 Yb8i_l1 = reshape(Yb8i_l1,V_segR.dim);
                %                 Yb8i_l2 = reshape(Yb8i_l2,V_segR.dim);
                %                 if Pwl==1
                %                     norm_ms = (1-Svx*sum_Jout)/(Vvx*sum(sum(sum(msP.*Yb8i_l1))));
                %                 else
                %                     norm_ms = (1-Svx*sum_Jout)/(Vvx*sum(sum(sum(msP.*Yb8i_l2))));
                %                 end
                %                 msP = msP*norm_ms;
                
            case 2 % tMCimg
                fid=fopen(f{fi,:},'rb');
                ms = fread(fid,V_segR.dim(1)*V_segR.dim(2)*V_segR.dim(3),'float32');
                fclose(fid);
                % Bonnery from Boas et al.
%                 sum_Jout = sum(ms(ms<0)/cs.par.nphotons);
%                 Svx = 5*5;
%                 Vvx = 5*5*5;
%                 [indices,~] = find(ms>0);
%                 msP = zeros(size(ms));
%                 msP(indices,1)=ms(indices,1);
%                 if Pwl==1
%                     norm_ms = (1-Svx*sum_Jout)/(Vvx*sum(msP.*Yb8i_l1));
%                 else
%                     norm_ms = (1-Svx*sum_Jout)/(Vvx*sum(msP.*Yb8i_l2));
%                 end
%                 msP = msP*norm_ms;
                %%%% version 2
               % ms = ms/cs.par.nphotons;
                sum_Jout = sum(ms(ms<0)/cs.par.nphotons);
                Svx = 5*5;
                Vvx = 5*5*5;
                [indices,~] = find(ms>0);
                msP = zeros(size(ms));
                msP(indices,1)=ms(indices,1);
                if Pwl==1
                    norm_ms = (1-Svx*abs(sum_Jout))/(Vvx*sum(msP.*Yb8i_l1));
                else
                    norm_ms = (1-Svx*abs(sum_Jout))/(Vvx*sum(msP.*Yb8i_l2));
                end
                msP = msP*norm_ms;
        end
        
        D_Sn =  unique(Cid(3,Cid(2,:)== Sn));% Detectors seeing source Sn
        for i = 1:size(D_Sn,2)% overview of detectors seen by source thanks to Cid
            for j = 1:size(fD,1)% only if detector has been selected by user...
                [~,fD_n,~] = fileparts(fD{j,:});
                
                if D_Sn(1,i) < 10, Ds_Sn = ['0' int2str(D_Sn(1,i))]; else Ds_Sn = int2str(D_Sn(1,i));end
                Dfn = ['D_No' Ds_Sn '_' int2str(wl(Pwl)) 'nm'];
                b = strcmp(fD_n,Dfn);
                if b % one pair found : processing goes on !
                    c = find(Cid(2,:)== Sn & Cid(3,:)==D_Sn(i) & Cwl==Pwl);
                    disp(['Sn=' int2str(Sn) ' et D_Sn=[' int2str(D_Sn) '] et wl=' int2str(Pwl) ' , DONC : ' int2str(c)]);
                    
                    switch cs.alg
                        case 1 % MCX
                            md=loadmc2(fullfile(cs_dir,[Dfn Oe]),V_segR.dim);
                            % Bonnery from Boas et al.
                            %                             sum_Jout = sum(md(md<0)/cs.par.nphotons);
                            %                             Svx = 5*5;
                            %                             Vvx = 5*5*5;
                            %                             mdP = md;
                            %                             mdP(mdP<0)=0;
                            %                             Yb8i_l1 = reshape(Yb8i_l1,V_segR.dim);
                            %                             Yb8i_l2 = reshape(Yb8i_l2,V_segR.dim);
                            %                             if Pwl==1
                            %                                 norm_md = (1-Svx*sum_Jout)/(Vvx*sum(sum(sum(mdP.*Yb8i_l1))));
                            %                             else
                            %                                 norm_md = (1-Svx*sum_Jout)/(Vvx*sum(sum(sum(mdP.*Yb8i_l2))));
                            %                             end
                            %                             mdP = mdP*norm_md;
                            
                            % sur la moyenne des points autour... : a coder
                            rayon = 3;
                            cntr = ceil(rayon/2);
                            phi0_S =0;
                            phi0_D =0;
                            for ir = 1:rayon
                                phi0_S = phi0_S + ms(cs.Pfp_rmiv(1,Sn)+(ir-cntr),cs.Pfp_rmiv(2,Sn)+(ir-cntr),cs.Pfp_rmiv(3,Sn)+(ir-cntr));
                                phi0_D = phi0_D + md(cs.Pfp_rmiv(1,D_Sn(ir))+(ir-cntr),cs.Pfp_rmiv(2,D_Sn(ir))+(ir-cntr),cs.Pfp_rmiv(3,D_Sn(ir))+(ir-cntr));
                            end
                            
                            % Sensitivity matrix
                            sens_sd = ms.*md / (1+(phi0_S + phi0_D)/2);
                            sens(isd+1,:) = reshape(sens_sd,[numel(sens_sd),1]);
                            isd = isd+1;
                            C = [C c]; % channels in the sensitivity matrix
                            
                        case 2 % tMCimg
                            fid=fopen(fullfile(cs_dir,[Dfn Oe]),'rb');
                            md = fread(fid,V_segR.dim(1)*V_segR.dim(2)*V_segR.dim(3),'float32');
                            fclose(fid);
                            % Bonnery from Boas et al.
                            sum_Jout = sum(md(md<0)/cs.par.nphotons);
                            Svx = 5*5;
                            Vvx = 5*5*5;
                            [indices,~] = find(md>0);
                            mdP = zeros(size(md));
                            mdP(indices,1)=md(indices,1);
                            if Pwl==1
                                norm_md = (1-Svx*abs(sum_Jout))/(Vvx*sum(mdP.*Yb8i_l1));
                            else
                                norm_md = (1-Svx*abs(sum_Jout))/(Vvx*sum(mdP.*Yb8i_l2));
                            end
                            mdP = mdP*norm_md;
                            %
                            
                            msR = reshape(msP,V_segR.dim);
                            mdR = reshape(mdP,V_segR.dim);
                            % sur la moyenne des points autour...
                            rayon = 1;
                            cntr = ceil(rayon/2);
                            phi0_S =0;
                            phi0_D =0;
                            for i = 1:rayon
                                phi0_S = phi0_S + msR(cs.Pfp_rmiv(1,Sn)+(i-cntr),cs.Pfp_rmiv(2,Sn)+(i-cntr),cs.Pfp_rmiv(3,Sn)+(i-cntr));
                                phi0_D = phi0_D + mdR(cs.Pfp_rmiv(1,D_Sn(i))+(i-cntr),cs.Pfp_rmiv(2,D_Sn(i))+(i-cntr),cs.Pfp_rmiv(3,D_Sn(i))+(i-cntr));
                            end
                            % Sensitivity matrix
                            sens_sd = msP.*mdP / (1+(phi0_S + phi0_D)/2);
                            sens(isd+1,:) = sens_sd';
                            isd = isd+1;
                            C = [C c]; % channels in the sensitivity matrix
                    end
                end
            end
        end
    end
end
NIRS.Cs.mcs{ics}.C = C;
save(job.NIRSmat{1,1},'NIRS');

sens_index = [C' sens];
sensC = sortrows(sens_index);
sensC(isnan(sensC))=0;
disp([int2str(length(isnan(sensC))) ' NaNs have been corrected in sensitivity matrix (' int2str(numel(sensC)) ' values)']);
sens = sensC(:,2:end);
save(fullfile(cs_dir,'sens.mat'),'sens');

sens_reshaped = zeros(V_segR.dim);
for i=1:size(sens,1)
    sens_reshaped = sens_reshaped + reshape(sens(i,:),V_segR.dim);
end
save(fullfile(cs_dir,'sensReshaped.mat'),'sens_reshaped');

%%% a mettre en option
V = struct('fname',fullfile(cs_dir,['sens' '.nii']),...
    'dim',  V_segR.dim,...
    'dt',   V_segR.dt,...
    'pinfo',V_segR.pinfo,...
    'mat',  V_segR.mat);

V = spm_create_vol(V);
spm_write_vol(V,log(sens_reshaped));
%%%

m_c1 = zeros(size(Yb8i));
m_c1(Yb8i==1)=1;

PVE = zeros(size(sens,1),2);
sens_sd_c1i = zeros(1,size(sens,2));
for i=1:size(sens,1)
    sens_sd = reshape(sens(i,:),V_segR.dim);
    V = struct('fname',fullfile(cs_dir,['banane_' int2str(sensC(i,1)) '.nii']),...
        'dim',  V_segR.dim,...
        'dt',   V_segR.dt,...
        'pinfo',V_segR.pinfo,...
        'mat',  V_segR.mat);
    
    V = spm_create_vol(V);
    spm_write_vol(V,sens_sd);
    
    % PVE : partial volume effect, ratio of optical densities in a whole
    % head compared to in the gray matter for the same activation
    for j=1:size(sens,2)
        if m_c1(j,1)==1
            sens_sd_c1i(1,j) = sens(i,j);
        end
    end
    PVE(i,1) = sensC(i,1);
    PVE(i,2) = sum(sens(i,:))/sum(sens_sd_c1i);%length(find(sens_sd_c1i>200))/length(find(sens(i,:)>200));
    
    V_c1 = struct('fname',fullfile(cs_dir,['banane_c1_' int2str(sensC(i,1)) '.nii']),...
        'dim',  V_segR.dim,...
        'dt',   V_segR.dt,...
        'pinfo',V_segR.pinfo,...
        'mat',  V_segR.mat);
    
    V_c1 = spm_create_vol(V_c1);
    spm_write_vol(V_c1,reshape(sens_sd_c1i,V_segR.dim));
    
end
save(fullfile(cs_dir,'PVE.mat'),'PVE');
out = 1;
end