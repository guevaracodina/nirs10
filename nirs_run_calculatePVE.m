function out = nirs_run_calculatePVE(job)
%%%%%%%%
% les Inf sont quand la lumiere ne passe pas par les zones du masque (BOLD par exemple)
%%%%%%%%
% Calcul du facteur de volume partiel
%
%_________________________________________________________________________
% Clément Bonnéry June 2011
% Methode from Hiraoka et al. 1993, Optical pathlength in inhomogeneous
% tissue

% valeurs donnees par Duncan :
% DPF 690 = 5.38 + 0.049 A 0.877
% DPF 744 = 5.11 + 0.106 A 0.723
% DPF 807 = 4.99 + 0.067 A 0.814
% DPF 832 = 4.67 + 0.062 A 0.819
DelPreviousData  = job.DelPreviousData;
try
    NewNIRSdir = job.NewDirCopyNIRS.CreateNIRSCopy.NewNIRSdir;
    NewDirCopyNIRS = 1;
catch
    NewDirCopyNIRS = 0;
end
% Loop over subjects
for iSubj=1:size(job.NIRSmat,1)
    % Load NIRS.mat
    try
        NIRS = [];
        load(job.NIRSmat{iSubj,1});
        
        NIRS.Dt.fir.DPF =[];
        NIRS.Dt.fir.PVE =[];
        
        % Get "cs" (current simulation) info
        if ~isempty(job.dir_in{1,1})
            cs_dir =  fileparts(job.dir_in{1,1});
            [dummy cs_ldir] = fileparts(cs_dir);
            ics =1;
            while ~strcmp(cs_ldir,NIRS.Cs.n{ics})
                ics =ics+1;
            end
        else
            [cs_dir dummy dummy1] = fileparts(job.NIRSmat{iSubj,1});
            ics = length(NIRS.Cs.n); % Use last simulation run
        end
        cs = NIRS.Cs.mcs{ics};
        if cs.alg==1, Oe='.mch'; elseif cs.alg==2, Oe='.his';end
        
        % Select history/mch files
        [fchar,dummy] = spm_select('FPList',cs_dir,Oe);
        for i=1:size(fchar,1)
            f{i,1} = deblank(fchar(i,:));
        end
        % Opt_ppts% Get optical properties of each layer of the medium and of the perturbatio
        if isfield(cs,'mu_subj') && isfield(cs.mu_subj,'muTRS')
            [dummy namef] = fileparts(f{1,1});
            num = str2num(namef(5:6));
            outOP = GetOpt_ppts('wl','D:\Users\Clément\DPF_testDuncan3\Opt_ppts_44sujets_MichP2S.mat',num);
        else
            outOP = GetOpt_ppts('wl');
        end
        opt_ppts = outOP{1};
        opt_ppts_perturb = outOP{2};
        
        % Read simulation outputs from his/mch file
        switch cs.alg
            case 1 %MCX
                count=0;
                for i=1:size(f,1)
                    [history, header]=loadmch(f{i,1});
                    % [history, header]=loadmch([f{i,1}(1:end-1) 'h']);
                    %                     [tk_dir, tk_file, tk_ext] = fileparts([f{i,1}(1:end-1) 'h']);
                    [tk_dir, tk_file, tk_ext] = fileparts(f{i,1});
                    %                     tk_wl = str2double(tk_file(end-4:end-2));
                    count = count+1;
                    History{count,1} = tk_file; % name of file
                    History{count,3} = header; % header info
                    History{count,2} = history; % for each detected photon:
                    % det # / number of scattering events / pathlength
                    % through each layer (nlayers+2 columns)
                end
                
            case 2 %tMCimg
                if strcmp(job.dir_in{:}(end),filesep)
                    job.dir_in{:} = job.dir_in{:}(1:end-1);
                end
                [t,dummy] = spm_select('FPList',job.dir_in{:},'.his');
                count =0;
                
                for k1=1:size(t,1)
                    tk =t(k1,:);
                    [tk_dir, tk_file, tk_ext] = fileparts(tk);
                    
                    if strcmp(tk_file(1),'S')%seulement si on a affaire a une source
                        tk_S = str2num(tk_file(5:6));
                        tk_wl = str2num(tk_file(end-4:end-2));
                        
                        code='mich';
                        switch code
                            case 'boas'
                                Cphore(1).Conc = 'HbO';
                                Cphore(1).Conc =60e-6;
                                Cphore(2).Name = 'Hb';
                                Cphore(2).Conc =40e-6;
                                [mua] = GetMua(tk_wl,Cphore);
                                
                                [fluence, nPhoton] = lirehis_clm(tk_file,muas,1);
                                
                            case 'mich'
                                ntissues = 6;
                                [history] = lirehis_clm(t(k1,:), cs.par.nphotons,cs.NDkpt,ntissues,cs.par.numTimeGates);
                                count = count+1;
                        end
                        History{count,1} = tk_file;
                        History{count,2} = history;
                    end
                end
            otherwise
                disp('The algorithm with which the simulation has been runned is not recognised.')
        end
        clear tk_wl tk_file tk_dir tk_ext history header
        
        Cid = NIRS.Cf.H.C.id;
        % 3 x NC : channel ID; source #; det #
        
        %Cgp (source-detector distance for each channel) : measure of the distance on the head
        Cgp = NIRS.Cf.H.C.gp;

        b0=0;b1=0;
        % b0 : nombre de photons comptés sur le détecteur en cours
        % b1 : nombre total de photons comptés
        
        % For each channel
        for Ci=1:NIRS.Cf.H.C.N
            try
                S_Ci = unique(Cid(2,Cid(1,:)== Ci));
                D_Ci = unique(Cid(3,Cid(1,:)== Ci));
                Cbloup = sum( (1:length(Cid(1,:))) .* (Cid(1,:) == Ci) );% vrai numero du canal
                % SOURCE
                if S_Ci < 10, S_Cin = ['0' int2str(S_Ci)]; else S_Cin = int2str(S_Ci);end
                tk_wl = NIRS.Cf.dev.wl(NIRS.Cf.H.C.wl(Cbloup)); % wavelength for this channel
                tk_n = ['S_No' S_Cin '_' int2str(tk_wl) 'nm'];
                
                iH=1;
                while strcmp(tk_n,History{iH,1})==0
                    iH=iH+1;
                end
                tk_Ci = iH;% indice de la bonne ligne pour le canal Cbloup
                
                % maintenant les detecteurs sont ranges comme dans le fichier de config, il y en a le
                % nombre de P-1 et cest dans lordre en enlevant la bonne source
                %%% on veut D_Ci
                %%%% on cherche le nombre total de sources et on enleve
                %%%% 1 :
                offset_cfgfile = cs.NSkpt-1;
                D_Ci_csPkpt = sum((1:size(cs.Pkpt,2)).*(cs.Pkpt==D_Ci+cs.NSinit));
                D_Ci_cfg = offset_cfgfile+D_Ci_csPkpt-cs.NSkpt;
                
                idx = History{tk_Ci,2}(:,1)==D_Ci_cfg;
                % 1 ligne par photon; idx recense les lignes pour lesquelles c'est
                %  le détecteur D_Ci_cfg qui a compté le photon
                
                % DETECTOR
                try
                    if D_Ci < 10, D_Cin = ['0' int2str(D_Ci)]; else D_Cin = int2str(D_Ci);end
                    tkD_n = ['D_No' D_Cin '_' int2str(tk_wl) 'nm'];
                    
                    iH=1;
                    while strcmp(tkD_n,History{iH,1})==0
                        iH=iH+1;
                    end
                    tkD_Ci = iH;% indice de la bonne ligne pour le canal Ci
                    
                    % sources are sorted the same way in the config file and in
                    % the NIRS matrix (in Pkpt) as they are sorted first
                    S_Ci_csPkpt = sum((1:size(cs.Pkpt,2)).*(cs.Pkpt==S_Ci));
                    S_Ci_cfg = S_Ci_csPkpt;
                    
                    idxD = History{tkD_Ci,2}(:,1)== S_Ci_cfg;
                    % 1 ligne par photon; idx recense les lignes pour lesquelles c'est
                    %  la source S_Ci_cfg qui a compté le photon
                catch
                    idxD = 0;
                end

                if sum(idx)==0 && sum(idxD)==0
                    PVF(1,Ci)=0;
                    DPF(1,Ci)=0;
                else
                    % on prepare les calculs puisque des photons ont ete
                    % detectes par le canal
                    if tk_wl==690
                        %                         muas = [opt_ppts{1,1}(1,1),opt_ppts{1,2}(1,1),opt_ppts{1,3}(1,1),opt_ppts{1,4}(1,1),opt_ppts{1,5}(1,1),opt_ppts_perturb{1}(1,1)];
                        muas = [opt_ppts(1,1,1),opt_ppts(1,2,1),opt_ppts(1,3,1),opt_ppts(1,4,1),opt_ppts(1,5,1),opt_ppts(1,6,1)];
                    elseif tk_wl ==830
                        %                         muas = [opt_ppts{2,1}(1,1),opt_ppts{2,2}(1,1),opt_ppts{2,3}(1,1),opt_ppts{2,4}(1,1),opt_ppts{2,5}(1,1),opt_ppts_perturb{2}(1,1)];
                        muas = [opt_ppts(2,1,1),opt_ppts(2,2,1),opt_ppts(2,3,1),opt_ppts(2,4,1),opt_ppts(2,5,1),opt_ppts(2,6,1)];
                    end
                    
                    % Source-detector distance
                    chord = norm(cs.P.Pfp_rmm(:,D_Ci_csPkpt)-cs.P.Pfp_rmm(:,S_Ci),2);
                    
                    test=1;
                    switch test
                        case 1
                            W0 = 1; % Okada : 1 for and to be calculated otherwise but it disappears anyway
                            if cs.nummed==6 % 6 layers (no perturbation)
                                L_Vi_Vphts = History{tk_Ci,2}(idx,3:size(History{tk_Ci,2},2));
                                L_Vi_Vphts = [L_Vi_Vphts; History{tkD_Ci,2}(idxD,3:size(History{tkD_Ci,2},2))];
                                %                                 elseif cs.nummed==11
                                %                                 L_Vi_Vphts = History{tk_Ci,2}(idx,3:8)+[History{tk_Ci,2}(idx,9:13),zeros(size(History{tk_Ci,2}(idx,9:13),1),1)];
                            elseif cs.nummed==12 % 12 layers (medium 6 layers + perturbation)
                                L_Vi_Vphts = History{tk_Ci,2}(idx,3:8)+History{tk_Ci,2}(idx,9:14);
                                L_Vi_Vphts = [L_Vi_Vphts; History{tkD_Ci,2}(idxD,3:8)+History{tkD_Ci,2}(idx,9:14)];
                            end
                            % muas' : nLayer x 1
                            % L_Vi_Vphts : nPhotons counted by this det x nLayer
                            W_phts = W0*exp(-sum(L_Vi_Vphts*muas',2));% one line per photon counted
                            b0 = b0+size(L_Vi_Vphts,1); % nPhotons counted by this det
                            b1 = b1+size(History{tk_Ci,2},1)+size(History{tkD_Ci,2},1); % total nPhotons detected
                            numerateur = zeros(1,6);
                            for ilayer=1:6
                                numerateur(1,ilayer) = sum(L_Vi_Vphts(:,ilayer).*W_phts,1);% sum over the photons
                            end
                            
                            PDP_Vi = numerateur/sum(W_phts);
                            % Differential pathlength factor,
                            % dimensionless, is mean pathglength /
                            % source-detector distance
                            PDPF_Vi(1:6,Ci) = PDP_Vi/chord;
                            
                            if cs.nummed==12
                                % Parcours moyen dans la zone BOLD
                                L_disturb_Vphts = History{tk_Ci,2}(idx,9:14);
                                L_disturb_Vphts = [L_disturb_Vphts; History{tkD_Ci,2}(idxD,9:14)];
                                
                                W_phts = W0*exp(-sum(L_disturb_Vphts*muas',2));
                                for ilayer=1:6
                                    num_disturb(1,ilayer) = sum(L_disturb_Vphts(:,ilayer).*W_phts,1);% sommation des photons
                                end
                                PDP_disturb = num_disturb/sum(W_phts);
                                PDPF_disturb(1:6,Ci) = PDP_disturb/chord;
                            end
                    end
                    disp(['DPF and PVE calculated for channel : ' int2str(Ci)]);
                end
            catch % le fichier d histoire n a pas ete trouve...
                disp(['DPF and PVE failed for channel (must probable reason : source is not part of the ROI): ' int2str(Ci)]);
                PDPF_Vi(:,Ci) =-10;
            end
        end
        if cs.nummed==12
            % PVF partial volume factor (PVF := DPF/PPF);
            PVF = sum(PDPF_Vi,1)./sum(PDPF_disturb,1);
            PVF(isinf(PVF)) = 0;
            save(fullfile(cs_dir,'PVF.mat'),'PVF');
            NIRS.Dt.fir.PVF = PVF;
        end
        
        disp(['b0 : ' int2str(b0) ' contre un total de ' int2str(b1)]);
        save(fullfile(cs_dir,'PDPF.mat'),'PDPF_Vi');
        NIRS.Dt.fir.PDPF = PDPF_Vi;
        
        % ----------------------------------------------------------------------- %
        if NewDirCopyNIRS
            [dirN fil1 ext1] =fileparts(job.NIRSmat{iSubj,1});
            dir2 = [dirN filesep NewNIRSdir];
            if ~exist(dir2,'dir'), mkdir(dir2); end;
            newNIRSlocation = fullfile(dir2,'NIRS.mat');
            save(newNIRSlocation,'NIRS');
            job.NIRSmat{iSubj,1} = newNIRSlocation;
        else
            save(job.NIRSmat{iSubj,1},'NIRS');
        end
        
        %         if DelPreviousData
        %             delete(rDtp{f,1});
        %         end
        
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        disp(['Calculus of PVE failed for subject ' int2str(iSubj)]);
    end
end
out.NIRSmat = job.NIRSmat;
end