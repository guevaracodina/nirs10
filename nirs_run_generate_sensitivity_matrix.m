function out = nirs_run_generate_sensitivity_matrix(job)
% Generates sensitivity matrix from output MC simulations files (MCX or
% tMCimg)
% FORMAT nirs_run_generate_sensitivity_matrix(outMCfiles,NIRSmat)
% outMCfiles - MC simulations output files (2pt or mc2)
% NIRSmat    - NIRS matrix
%_______________________________________________________________________
%
% First the MC output is normalised according to Boas et al. and Fang et
% al.
% As some simulations may have failed the code is trying to recover all the
% pairs. One loop is applied to sources then the paired detectors are
% searched.
% After normalisation, PVE is calculated.
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire
% Clément Bonnéry 03/2011

for Idx=1:size(job.NIRSmat,1)
    %Load NIRS.mat information
    try
        clear sens
        [NIRS newNIRSlocation]= nirs_load(job.NIRSmat{Idx,1},job.NIRSmatCopyChoice,job.force_redo);
        job.NIRSmat{Idx,1} = newNIRSlocation;
        if ~isempty(NIRS) && (~isfield(NIRS.flags,'sensOK') || job.force_redo)
            if ~isempty(job.outMCfiles{1,1})
                f = job.outMCfiles;
                cs_dir =  fileparts(f{1,:});
                [dummy cs_ldir] = fileparts(cs_dir);
                %             cs_ldir = cs_dir(max(strfind(cs_dir,'\'))+(length('MC')+1):end);
                
                ics =1;
                while ~strcmp(cs_ldir,NIRS.Cs.n{ics})
                    ics =ics+1;
                end
                cs = NIRS.Cs.mcs{ics};
                if cs.alg==1, Oe='.mc2'; elseif cs.alg==2, Oe='.2pt';end
                
            else
                [cs_dir dummy dummy1] = fileparts(job.NIRSmat{Idx,1});
                ics = length(NIRS.Cs.n);
                cs_ldir = NIRS.Cs.n{1,ics};
                cs = NIRS.Cs.mcs{ics};
                if cs.alg==1, Oe='.mc2'; elseif cs.alg==2, Oe='.2pt';end
                
                [fchar,dummy] = spm_select('FPList',cs_dir,Oe);
                for i=1:size(fchar,1)
                    f{i,1} = fchar(i,:);
                end
            end
            
            %%%Opt_ppts
            % % % % % % % outOP = GetOpt_ppts('wl');
            if isfield(cs,'mu_subj') && isfield(cs.mu_subj,'muTRS')
                [dummy namef] = fileparts(f{1,1});
                num = str2num(namef(5:6));
                outOP = GetOpt_ppts('wl','D:\Users\Clément\DPF_testDuncan3\Opt_ppts_44sujets_MichP2S.mat',num);
            else
                outOP = GetOpt_ppts('wl');
            end
            opt_ppts = outOP{1};
            opt_ppts_perturb = outOP{2};
            
            try MCX_g = NIRS.Cs.mcs{ics}.MCX_g; catch, MCX_g =1; end
            
            
            segR = cs.segR;
            V_segR = spm_vol(segR);
            
            % lecture du fichier binaire pour avoir les mu_a
            b8i = cs.b8i;
            fid=fopen(b8i,'rb');
            Yb8i = fread(fid,prod(V_segR.dim),'uint8');
            fclose(fid);
            % : mua mus' ... on prend les mu_a
            Yb8i_l1 = zeros(size(Yb8i));
            Yb8i_l1(Yb8i==1)= opt_ppts(1,1,1);
            Yb8i_l1(Yb8i==2)= opt_ppts(1,2,1);
            Yb8i_l1(Yb8i==3)= opt_ppts(1,3,1);
            Yb8i_l1(Yb8i==4)= opt_ppts(1,4,1);
            Yb8i_l1(Yb8i==5)= opt_ppts(1,5,1);
            Yb8i_l1(Yb8i==6)= opt_ppts(1,6,1);
            
            Yb8i_l2 = zeros(size(Yb8i));
            Yb8i_l2(Yb8i==1)= opt_ppts(2,1,1);
            Yb8i_l2(Yb8i==2)= opt_ppts(2,2,1);
            Yb8i_l2(Yb8i==3)= opt_ppts(2,3,1);
            Yb8i_l2(Yb8i==4)= opt_ppts(2,4,1);
            Yb8i_l2(Yb8i==5)= opt_ppts(2,5,1);
            Yb8i_l2(Yb8i==6)= opt_ppts(2,6,1);
            
            if size(Yb8i_l1,1)==1%cs.alg==2
                Yb8i_l1 = Yb8i_l1';
                Yb8i_l2 = Yb8i_l2';
            end
            
            wl = NIRS.Cf.dev.wl; %= [830 690];
            Cid = NIRS.Cf.H.C.id;
            Cwl = NIRS.Cf.H.C.wl;
            C =[];
            NfD = 0;
            NfS = 0;
            isd =0; % count of lignes in sensitivity matrix
            
            Skpt = cs.Pkpt(1:cs.NSkpt);
            Dkpt = cs.Pkpt(cs.NSkpt+1:cs.NSkpt+cs.NDkpt)-8;
            
            %%%%%%%%%%
            for fi = 1:size(f,1)
                [dummy,fn,dummy2] = fileparts(f{fi,:});
                
                if strcmp(fn(1:1),'D')% stores all ''D_No.._...nm.2pt'' files in a list
                    fD{NfD+1,:} = f{fi,:};
                    NfD = NfD+1;
                    
                elseif strcmp(fn(1:1),'S')% same with ''S_No.._...nm.2pt'' files
                    fS{NfS+1,:} = f{fi,:};
                    NfS = NfS+1;
                end
            end
            %%%% bug dans le cas ou on a une longueur d onde qui a buggee
            if rem(NfD,2)
                %%% on cherche quelle sim a echouee et on supprime le point +
                %%% message
                disp('a faire')
            else
                ND = NfD/2;
            end
            if rem(NfS,2)
                %%% on cherche quelle sim a echouee et on supprime le point +
                %%% message
                disp('a faire')
            else
                NS = NfS/2;
            end
            
            %Pre-define the size of the sensitivity matrix for the memory
            sens = zeros(length(Cwl),length(Yb8i));
            
            %added to show progress
            disp('sens matrix pre-created');
            
            for fiS = 1:NfS
                [dummy,fSn,dummy2] = fileparts(fS{fiS,:});
                [dummy,Pwl] = find(wl==str2double(fSn(8:10)));
                Sn = str2double(fSn(5:6));
                switch cs.alg
                    case 1 % MCX
                        % already normalized (http://mcx.sourceforge.net/cgi-bin/index.cgi?Doc/README : 6.1 output files)
                        try %%%%%%%PROVISOIRE PROVISOIRE PROVISOIRE
                            numgates = min(cs.numTimeGates,MCX_g);
                        catch
                            numgates = MCX_g;
                        end
                        
                        ms4=loadmc2(fS{fiS,:},[V_segR.dim numgates],'float');
                        if numgates~=1
                            ms=sum(ms4,4);
                        else
                            ms=ms4;
                        end
                        
                    case 2 % tMCimg
                        fid=fopen(fS{fiS,:},'rb');
                        %                     T = spm_type(fid,'bits');
                        ms = fread(fid,V_segR.dim(1)*V_segR.dim(2)*V_segR.dim(3),'float64');%,'float32');
                        fclose(fid);
                        %%%%%%%%%%%%%%%%%%%%%%%
                        % % % % %                     code ;ich sur les fluences voir aussi le code
                        % de Boas
                        % % % % %                                                     % Fluence au détecteur selon éq. (3) de Boas 2002
                        % % % % %                                 fluence(pi) = sum(photons_weight)./nbPhotonsTotSimu;
                        %%%%%%%%%%%%%%%%%%%%%%%
                        sum_Jout = sum(ms(ms<0)/(cs.par.nphotons));
                        Svx = 5*5;
                        Vvx = 5*5*5;
                        [indices,~] = find(ms>0);
                        msP = zeros(size(ms));
                        msP(indices,1)=ms(indices,1);
                        if Pwl==1
                            norm_ms = (1-Svx*sum_Jout)/(Vvx*sum(msP.*Yb8i_l1));
                        else
                            norm_ms = (1-Svx*sum_Jout)/(Vvx*sum(msP.*Yb8i_l2));
                        end
                        msP = msP*norm_ms;
                end
                
                D_Sn =  unique(Cid(3,Cid(2,:)== Sn));% Detectors seeing source Sn
                %Assign array
                %sens = zeros(NfS,length(Yb8i));
                %sens = zeros(length(Cwl),length(Yb8i));
                for i = 1:size(D_Sn,2)% overview of detectors seen by source according to Cid
                    for j = 1:NfD% only if detector has been selected by user...
                        [dummy,fD_n,dummy2] = fileparts(fD{j,:});
                        
                        if D_Sn(1,i) < 10, Ds_Sn = ['0' int2str(D_Sn(1,i))]; else Ds_Sn = int2str(D_Sn(1,i));end
                        Dfn = ['D_No' Ds_Sn '_' int2str(wl(Pwl)) 'nm'];
                        b = strcmp(fD_n,Dfn);
                        if b % one pair found : processing goes on !
                            c = find(Cid(2,:)== Sn & Cid(3,:)==D_Sn(i) & Cwl==Pwl);
                            %disp(['Sn=' int2str(Sn) ' et D_Sn=[' int2str(D_Sn) '] et wl=' int2str(Pwl) ' , DONC : ' int2str(c)]);
                            
                            D_Pktp = NS+sum((1:length(Dkpt)).*(Dkpt==D_Sn(i)));
                            S_Pkpt = sum((1:length(Skpt)).*(Skpt==Sn));
                            switch cs.alg
                                case 1 % MCX
                                    % already normalized (http://mcx.sourceforge.net/cgi-bin/index.cgi?Doc/README : 6.1 output files)
                                    md4=loadmc2(fullfile(cs_dir,[Dfn Oe]),[V_segR.dim numgates],'float');
                                    if numgates~=1
                                        md=sum(md4,4);
                                    else
                                        md=md4;
                                    end
                                    
                                    %%% calcul de phi0
                                    % boule autour de la position du point P
                                    vxr = 3;%voxel radius
                                    
                                    ms_N = ms(round(max(cs.Pfp_rmiv(1,S_Pkpt)-vxr,1)):round(min(cs.Pfp_rmiv(1,S_Pkpt)+vxr,size(ms,1))),...
                                        round(max(cs.Pfp_rmiv(2,S_Pkpt)-vxr,1)):round(min(cs.Pfp_rmiv(2,S_Pkpt)+vxr,size(ms,2))),...
                                        round(max(cs.Pfp_rmiv(3,S_Pkpt)-vxr,1)):round(min(cs.Pfp_rmiv(3,S_Pkpt)+vxr,size(ms,3))));
                                    md_N = md(round(max(cs.Pfp_rmiv(1,D_Pktp)-vxr,1)):round(min(cs.Pfp_rmiv(1,D_Pktp)+vxr,size(md,1))),...
                                        round(max(cs.Pfp_rmiv(2,D_Pktp)-vxr,1)):round(min(cs.Pfp_rmiv(2,D_Pktp)+vxr,size(md,2))),...
                                        round(max(cs.Pfp_rmiv(3,D_Pktp)-vxr,1)):round(min(cs.Pfp_rmiv(3,D_Pktp)+vxr,size(md,3))));
                                    
                                    phi0_S = max(ms_N(:)); %%???
                                    phi0_D = max(md_N(:));
                                    
                                    % Sensitivity matrix
                                    sens_sd = ms.*md / ((phi0_S + phi0_D)/2);
                                    sens(isd+1,:) = reshape(sens_sd,[numel(sens_sd),1]);
                                    isd = isd+1;
                                    %add to show which step
                                    disp(['sens matrix treated for the colomn' num2str(isd)]);
                                    
                                    C = [C c]; % channels in the sensitivity matrix
                                    
                                case 2 % tMCimg
                                    fid=fopen(fullfile(cs_dir,[Dfn Oe]),'rb');
                                    md = fread(fid,V_segR.dim(1)*V_segR.dim(2)*V_segR.dim(3),'float64');%,'float32');
                                    fclose(fid);
                                    % Bonnery from Boas et al.
                                    sum_Jout = sum(md(md<0)/(cs.par.nphotons*10));
                                    Svx = 5*5; %???
                                    Vvx = 5*5*5;
                                    [indices,dummy] = find(md>0);
                                    mdP = zeros(size(md));
                                    mdP(indices,1)=md(indices,1);
                                    if Pwl==1
                                        norm_md = (1-Svx*sum_Jout)/(Vvx*sum(mdP.*Yb8i_l1));
                                    else
                                        norm_md = (1-Svx*sum_Jout)/(Vvx*sum(mdP.*Yb8i_l2));
                                    end
                                    mdP = mdP*norm_md;
                                    
                                    msR = reshape(msP,V_segR.dim);
                                    mdR = reshape(mdP,V_segR.dim);
                                    %%% calcul de phi0
                                    vxr = 3;%voxel radius in millimeters -- assumption to make more visible and explain
                                    
                                    ms_N = msR(max(cs.Pfp_rmiv(1,S_Pkpt)-vxr,1):min(cs.Pfp_rmiv(1,S_Pkpt)+vxr,size(msR,1)),...
                                        max(cs.Pfp_rmiv(2,S_Pkpt)-vxr,1):min(cs.Pfp_rmiv(2,S_Pkpt)+vxr,size(msR,2)),...
                                        max(cs.Pfp_rmiv(3,S_Pkpt)-vxr,1):min(cs.Pfp_rmiv(3,S_Pkpt)+vxr,size(msR,3)));
                                    
                                    md_N = mdR(max(cs.Pfp_rmiv(1,D_Pktp)-vxr,1):min(cs.Pfp_rmiv(1,D_Pktp)+vxr,size(mdR,1)),...
                                        max(cs.Pfp_rmiv(2,D_Pktp)-vxr,1):min(cs.Pfp_rmiv(2,D_Pktp)+vxr,size(mdR,2)),...
                                        max(cs.Pfp_rmiv(3,D_Pktp)-vxr,1):min(cs.Pfp_rmiv(3,D_Pktp)+vxr,size(mdR,3)));
                                    
                                    phi0_S = max(ms_N(:));
                                    phi0_D = max(md_N(:));
                                    
                                    % Sensitivity matrix
                                    sens_sd = msP.*mdP / ((phi0_S + phi0_D)/2);
                                    
                                    sens(isd+1,:) = sens_sd';
                                    isd = isd+1;
                                    C = [C c]; % channels in the sensitivity matrix
                            end
                        end
                    end
                end
            end
            
            %added to show which step
            disp('Before save sen matrix');
            
            %pre-define SensC to be prepared for memory
            %tempCsize = size(C);
            %sens_index = zeros(tempCsize(2), tempCsize(1)+ length(Yb8i));
            %sensC = zeros(tempCsize(2), tempCsize(1)+ length(Yb8i));
            
            NIRS.Cs.mcs{ics}.C = C;
            save(job.NIRSmat{1,1},'NIRS');
            
            sens_sparse = sparse(sens);
            clear sens;
            sens = sens_sparse;
            clear sens_sparse;
            
            [Csorted Cindex] = sortrows(C');
            sens = sens(Cindex,:);
            %sens_index = [C' sens];
            %sensC = sortrows(sens_index);
            %sensC(isnan(sensC))=0;
            %disp([int2str(sum(sum(sum(isnan(sensC))))) ' NaNs have been corrected in sensitivity matrix (' int2str(numel(sensC)) ' values)']);
            %sens = sensC(:,2:end);
            save_sens(cs_dir,sens);
            
            %save(fullfile(cs_dir,'sens.mat'),'sens','-v7.3'); %PP
            
            sens_reshaped = zeros(V_segR.dim);
            for i=1:size(sens,1)
                sens_reshaped = sens_reshaped + reshape(sens(i,:),V_segR.dim);
            end
            save(fullfile(cs_dir,'Sum_sensReshaped.mat'),'sens_reshaped','-v7.3'); %PP
            
            
            %added to show progress
            disp('after save sen matrix');
            
            %         sensR = reshape(sens,[size(sens,1) V_segR.dim]);
            %         for i=1:size(sens,1)
            %             sens2 = squeeze(sensR(i,:,:,:));
            %             %Save as a 4D .nii volume
            %             V(i) = nirs_create_vol(fullfile(cs_dir,['sens' int2str(i) '.nii']),...
            %                     V_segR.dim, [16,0], V_segR.pinfo, V_segR.mat, sens2);
            %         end
            %         %Save as a 4D .nii volume
            %         V4 = spm_file_merge(V,'sens.nii',0);
            %
            
            
            V2 = nirs_create_vol(fullfile(cs_dir,'log_sum_sens.nii'),...
                V_segR.dim, [16,0], V_segR.pinfo, V_segR.mat, log(sens_reshaped));
            
            %         m_c1 = zeros(size(Yb8i));
            %         m_c1(Yb8i==1)=1;
            %         NIRS.sens = V4;
            NIRS.flags.sensOK = 1;
            save(job.NIRSmat{Idx,1},'NIRS');
        end
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        disp(['Could not run MonteCarlo sensitivity matrix for subject' int2str(Idx) ' for ' job.NIRSmat{Idx,1}]);
    end
end
out.NIRSmat = job.NIRSmat;