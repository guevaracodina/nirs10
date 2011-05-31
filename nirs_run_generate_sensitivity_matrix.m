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
        NIRS = [];
        load(job.NIRSmat{Idx,1});

        if ~isempty(job.outMCfiles)
            f = job.outMCfiles;
            cs_dir =  fileparts(f{1,:});
            cs_ldir = cs_dir(max(strfind(cs_dir,'\'))+(length('MC')+1):end);
            
            ics =1;
            while ~strcmp(cs_ldir,NIRS.Cs.n{ics})
                ics =ics+1;
            end
        else
            ics = length(NIRS.Cs.n);
        end
        
        cs = NIRS.Cs.mcs{ics};
        try MCX_g = NIRS.Cs.mcs{ics}.MCX_g; catch, MCX_g =1; end
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
        NfD = 0;
        NfS = 0;
        isd =0; % count of lignes in sensitivity matrix

        %%%%%%%%%%
        for fi = 1:size(f,1)
            [~,fn,~] = fileparts(f{fi,:});

            if strcmp(fn(1:1),'D')% stores all ''D_No.._...nm.2pt'' files in a list
                fD{NfD+1,:} = f{fi,:};
                NfD = NfD+1;

            elseif strcmp(fn(1:1),'S')% same with ''S_No.._...nm.2pt'' files
                fS{NfS+1,:} = f{fi,:};
                NfS = NfS+1;
            end
        end
        ND = NfD/2;
        NS = NfS/2;

        for fiS = 1:NfS
            [~,fSn,~] = fileparts(fS{fiS,:});
            [~,Pwl] = find(wl==str2double(fSn(8:10)));
            Sn = str2double(fSn(5:6));
            switch cs.alg
                case 1 % MCX
                    % already normalized (http://mcx.sourceforge.net/cgi-bin/index.cgi?Doc/README : 6.1 output files)
                    numgates = min(cs.numTimeGates,MCX_g);

                    ms=loadmc2(fS{fiS,:},[V_segR.dim numgates],'float');
                    cw_mcx=sum(ms,4);
                    ms = cw_mcx;

                case 2 % tMCimg
                    fid=fopen(fS{fiS,:},'rb');
                    ms = fread(fid,V_segR.dim(1)*V_segR.dim(2)*V_segR.dim(3),'float32');
                    fclose(fid);
                    % Bonnery from Boas et al.
                    %                 sum_Jout = sum(ms(ms<0)/(cs.par.nphotons*10));
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
                    % ms = ms/(cs.par.nphotons*10);
                    sum_Jout = sum(ms(ms<0)/(cs.par.nphotons));
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
                for j = 1:NfD% only if detector has been selected by user...
                    [~,fD_n,~] = fileparts(fD{j,:});

                    if D_Sn(1,i) < 10, Ds_Sn = ['0' int2str(D_Sn(1,i))]; else Ds_Sn = int2str(D_Sn(1,i));end
                    Dfn = ['D_No' Ds_Sn '_' int2str(wl(Pwl)) 'nm'];
                    b = strcmp(fD_n,Dfn);
                    if b % one pair found : processing goes on !
                        c = find(Cid(2,:)== Sn & Cid(3,:)==D_Sn(i) & Cwl==Pwl);
                        disp(['Sn=' int2str(Sn) ' et D_Sn=[' int2str(D_Sn) '] et wl=' int2str(Pwl) ' , DONC : ' int2str(c)]);

                        switch cs.alg
                            case 1 % MCX
                                % already normalized (http://mcx.sourceforge.net/cgi-bin/index.cgi?Doc/README : 6.1 output files)
                                %                              md=loadmc2(fullfile(cs_dir,[Dfn Oe]),V_segR.dim,'float');
                                md=loadmc2(fullfile(cs_dir,[Dfn Oe]),[V_segR.dim numgates],'float');
                                cw_md=sum(md,4);
                                md = cw_mcx;

                                %%% calcul de phi0
                                % boule autour de la position du point P
                                vxr = 3;%voxel radius
                                ms_N = ms(cs.Pfp_rmiv(1,Sn)-vxr:cs.Pfp_rmiv(1,Sn)+vxr,cs.Pfp_rmiv(2,Sn)-vxr:cs.Pfp_rmiv(2,Sn)+vxr,cs.Pfp_rmiv(3,Sn)-vxr:cs.Pfp_rmiv(3,Sn)+vxr);
                                md_N = md(cs.Pfp_rmiv(1,NS+D_Sn(i))-vxr:cs.Pfp_rmiv(1,NS+D_Sn(i))+vxr,cs.Pfp_rmiv(2,NS+D_Sn(i))-vxr:cs.Pfp_rmiv(2,NS+D_Sn(i))+vxr,cs.Pfp_rmiv(3,NS+D_Sn(i))-vxr:cs.Pfp_rmiv(3,NS+D_Sn(i))+vxr);

                                phi0_S = max(ms_N(:));
                                phi0_D = max(md_N(:));
                                disp(['phi0=' int2str((phi0_S + phi0_D)/2)]);
                                % Sensitivity matrix
                                sens_sd = ms.*md / ((phi0_S + phi0_D)/2);
                                sens(isd+1,:) = reshape(sens_sd,[numel(sens_sd),1]);
                                isd = isd+1;
                                C = [C c]; % channels in the sensitivity matrix

                            case 2 % tMCimg
                                fid=fopen(fullfile(cs_dir,[Dfn Oe]),'rb');
                                md = fread(fid,V_segR.dim(1)*V_segR.dim(2)*V_segR.dim(3),'float32');
                                fclose(fid);
                                % Bonnery from Boas et al.
                                sum_Jout = sum(md(md<0)/(cs.par.nphotons*10));
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

                                msR = reshape(msP,V_segR.dim);
                                mdR = reshape(mdP,V_segR.dim);
                                %%% calcul de phi0
                                vxr = 3;%voxel radius
                                ms_N = ms(cs.Pfp_rmiv(1,Sn)-vxr:cs.Pfp_rmiv(1,Sn)+vxr,cs.Pfp_rmiv(2,Sn)-vxr:cs.Pfp_rmiv(2,Sn)+vxr,cs.Pfp_rmiv(3,Sn)-vxr:cs.Pfp_rmiv(3,Sn)+vxr);
                                md_N = ms(cs.Pfp_rmiv(1,NS+D_Sn(i))-vxr:cs.Pfp_rmiv(1,NS+D_Sn(i))+vxr,cs.Pfp_rmiv(2,NS+D_Sn(i))-vxr:cs.Pfp_rmiv(2,NS+D_Sn(i))+vxr,cs.Pfp_rmiv(3,NS+D_Sn(i))-vxr:cs.Pfp_rmiv(3,NS+D_Sn(i))+vxr);

                                phi0_S = max(ms_N(:));
                                phi0_D = max(md_N(:));
                                % Sensitivity matrix
                                phi0 = (phi0_S + phi0_D)/2;
                                disp(['phi0=' int2str(phi0)]);
                                sens_sd = msP.*mdP / phi0;
                                sens(isd+1,:) = sens_sd';
                                isd = isd+1;
                                C = [C c]; % channels in the sensitivity matrix
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
        disp([int2str(sum(isnan(sensC))) ' NaNs have been corrected in sensitivity matrix (' int2str(numel(sensC)) ' values)']);
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
            'dt',   [16,0],...
            'pinfo',V_segR.pinfo,...
            'mat',  V_segR.mat);

        V = spm_create_vol(V);
        spm_write_vol(V,log(sens_reshaped));
        V1 = struct('fname',fullfile(cs_dir,['ms1' '.nii']),...
            'dim',  V_segR.dim,...
            'dt',   [16,0],...
            'pinfo',V_segR.pinfo,...
            'mat',  V_segR.mat);
        V1 = spm_create_vol(V1);
        spm_write_vol(V1,log(ms));
        V2 = struct('fname',fullfile(cs_dir,['md1' '.nii']),...
            'dim',  V_segR.dim,...
            'dt',   [16,0],...
            'pinfo',V_segR.pinfo,...
            'mat',  V_segR.mat);
        V2 = spm_create_vol(V2);
        spm_write_vol(V2,log(md));
        %%%

        m_c1 = zeros(size(Yb8i));
        m_c1(Yb8i==1)=1;

        % sens_sd_c1i = zeros(1,size(sens,2));
        sens_c1 =  zeros(size(sens));
        sens_reshaped_c1 =  zeros(size(sens_reshaped));

        for i=1:size(sens,1)
            sens_sd = reshape(sens(i,:),V_segR.dim);
            V = struct('fname',fullfile(cs_dir,['banane_' int2str(sensC(i,1)) '.nii']),...
                'dim',  V_segR.dim,...
                'dt',   V_segR.dt,...
                'pinfo',V_segR.pinfo,...
                'mat',  V_segR.mat);

            V = spm_create_vol(V);
            spm_write_vol(V,sens_sd);

            for j=1:size(sens,2)
                if m_c1(j,1)==1
                    %             sens_sd_c1i(1,j) = sens(i,j);
                    sens_c1(i,j) = sens(i,j);
                end
            end
            sens_reshaped_c1 = sens_reshaped_c1 + reshape(sens_c1(i,:),V_segR.dim);

            V_c1 = struct('fname',fullfile(cs_dir,['banane_c1_' int2str(sensC(i,1)) '.nii']),...
                'dim',  V_segR.dim,...
                'dt',   V_segR.dt,...
                'pinfo',V_segR.pinfo,...
                'mat',  V_segR.mat);

            V_c1 = spm_create_vol(V_c1);
            spm_write_vol(V_c1,reshape(sens_c1(i,:),V_segR.dim));%sens_sd_c1i,V_segR.dim));

        end
        save(job.NIRSmat{Idx,1},'NIRS');       
    catch exception
        disp(exception.identifier);
        disp(['Could not run MonteCarlo sensitivity matrix for subject' int2str(Idx)]);
    end  
end
out.NIRSmat = job.NIRSmat;