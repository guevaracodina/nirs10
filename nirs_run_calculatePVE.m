function out = nirs_run_calculatePVE(job)
%%%%%%%%
% les Inf sont quand la lumiere ne passe pas par les zones du masque (BOLD par exemple)
%%%%%%%%
% Calcul du facteur de volume partiel
%
% BOlDdata;  (SEULEMENT SI ON VEUT EFFECTUER UN PROBLÈME DIRECT, i.e. PROJETER CES DONNÉES DANS L'ESPACE DES PAIRES)
% Masque du changement de signal BOLD (rBOLD) en % par rapport au baseline,
% matrice en voxels (interpolée à la résolution de l'IRM anatomique)
% size : nHbX*V x 1
% il faut que dConcfMRI
% soit une matrice (déroulée en un vecteur) de V éléments (pour chaque
% espèce HbX concaténée) qui donne un masque qui vaut deltaHbO et R dans
% tout le volume.
% De l'analyse Neurolens : FIR estimée normalisée * valeur moyenne
% sur la ROI * masque qui vaut 1 dans la ROI et 0 ailleurs, sommée (moyennée?) sur
% tmin:tmax.
%
% TRS : donne les données nécessaires au calcul du facteur de calibration du BOLD
% V0, OEF0 = (SO2in-SO2out)/SO2in = 1-ScO20, HbR0
%
% Quantités d'intérêt calculées :
% 	DPF (differential pathlength factor);
% 	PPF (partial pathlength factor dans la perturbation, si on a roulé une simulation
% avec perturbation (on suppose que le dernier tissu est une perturbation)***);
%
% V0 = 0.05; % venous blood volume fraction
% OEF0 = 1-0.6975; % baseline oxygen extraction fraction
% HbR0 = 131.*(1-0.758); % baseline HbR concentration (=HbT0*(1-SO20)) in uM
% TE = 30e-3; % 30 ms in our BOLD sequence
% nu0 = 80.6; % s^-1, at 3T
% alpha = -4.3 * nu0 * TE * V0 * OEF0 / HbR0; % Facteur de calibration du BOLD(relative change BOLD / uM) : BOLD(t)/BOLD(0) = alpha * [HbR(t)]
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


outOP = GetOpt_ppts('wl');
opt_ppts = outOP{1};
opt_ppts_perturb = outOP{2};

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
        
        
        if ~isempty(job.dir_in{1,1})
            f = job.outMCfiles;
            cs_dir =  fileparts(f{1,:});
            [dummy cs_ldir] = fileparts(cs_dir);
            
            ics =1;
            while ~strcmp(cs_ldir,NIRS.Cs.n{ics})
                ics =ics+1;
            end
            cs = NIRS.Cs.mcs{ics};
            if cs.alg==1, Oe='.mch'; elseif cs.alg==2, Oe='.his';end
            
        else
            [cs_dir dummy dummy1] = fileparts(job.NIRSmat{iSubj,1});
            ics = length(NIRS.Cs.n);
            cs_ldir = NIRS.Cs.n{1,ics};
            cs = NIRS.Cs.mcs{ics};
            if cs.alg==1, Oe='.mch'; elseif cs.alg==2, Oe='.his';end
            
            [fchar,dummy] = spm_select('FPList',cs_dir,Oe);
            for i=1:size(fchar,1)
                f{i,1} = fchar(i,:);
            end
        end
        
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
                    History{count,1} = tk_file;
                    History{count,3} = header;
                    History{count,2} = history;
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
        
        % cas caca, il faudra REFAIRE PLUS PROPRE
        Cid = NIRS.Cf.H.C.id;
        
        %Cgp
        Cgp = NIRS.Cf.H.C.gp;
        %mesure par les Pfp
        %Cgp =
        
        for Ci=1:NIRS.Cf.H.C.N
            try
                S_Ci = unique(Cid(2,Cid(1,:)== Ci));
                D_Ci = unique(Cid(3,Cid(1,:)== Ci));
                Cbloup = sum((1:length(Cid(1,:))).*(Cid(1,:) == Ci));
                %%% on reconstruit le nom
                if S_Ci < 10, S_Cin = ['0' int2str(S_Ci)]; else S_Cin = int2str(S_Ci);end
                tk_wl = NIRS.Cf.dev.wl(NIRS.Cf.H.C.wl(Cbloup));
                tk_n = ['S_No' S_Cin '_' int2str(tk_wl) 'nm'];
                %                 [dummy, tk_file, dummy1] = fileparts(tk_n);
                
                iH=1;
                while strcmp(tk_n,History{iH,1})==0
                    iH=iH+1;
                end
                tk_Ci = iH;% indice de la bonne ligne pour le canal Ci
                
                % maintenant les detecteurs sont ranges comme
                % dans le fichier de config, il y en a le
                % nombre de P-1 et cest dans lordre en
                % enlevant la bonne source
                %%% on veut D_Ci
                %%%% on cherche le nombre total de sources et on enleve
                %%%% 1 :
                offset_cfgfile = cs.NSkpt-1;
                D_Ci_csPkpt = sum((1:size(cs.Pkpt,2)).*(cs.Pkpt==D_Ci+cs.NSinit));
                D_Ci_cfg = offset_cfgfile+D_Ci_csPkpt-cs.NSkpt;
                
                idx = History{tk_Ci,2}(:,1)==D_Ci_cfg;
                if sum(idx)==0
                    PVF(1,Ci)=0;
                    DPF(1,Ci)=0;
                else
                    % on prepare les calculs puisque des photons ont ete
                    % detectes par le canal
                    if tk_wl==690
%                         muas = [opt_ppts{1,1}(1,1),opt_ppts{1,2}(1,1),opt_ppts{1,3}(1,1),opt_ppts{1,4}(1,1),opt_ppts{1,5}(1,1),opt_ppts_perturb{1}(1,1)];
                        muas = [opt_ppts{1,1}(1,1),opt_ppts{1,2}(1,1),opt_ppts{1,3}(1,1),opt_ppts{1,4}(1,1),opt_ppts{1,5}(1,1),opt_ppts{1,6}(1,1)];
                    elseif tk_wl ==830
%                         muas = [opt_ppts{2,1}(1,1),opt_ppts{2,2}(1,1),opt_ppts{2,3}(1,1),opt_ppts{2,4}(1,1),opt_ppts{2,5}(1,1),opt_ppts_perturb{2}(1,1)];
                        muas = [opt_ppts{2,1}(1,1),opt_ppts{2,2}(1,1),opt_ppts{2,3}(1,1),opt_ppts{2,4}(1,1),opt_ppts{2,5}(1,1),opt_ppts{2,6}(1,1)];
                    end
                    
                    chord = norm(cs.P.Pfp_rmm(:,D_Ci_csPkpt)-cs.P.Pfp_rmm(:,S_Ci),2);
                    
                    test=1;
                    switch test
                        case 1
                            % POUR MCX
                            %                     if b==0
                            W0 = 1;
                            %                     else
                            %                     end
                            if cs.nummed==6
                                L_Vi_Vphts = History{tk_Ci,2}(idx,3:size(History{tk_Ci,2},2));
%                                 elseif cs.nummed==12
%                                 L_Vi_Vphts = History{tk_Ci,2}(idx,3:8)+[History{tk_Ci,2}(idx,9:13),zeros(size(History{tk_Ci,2}(idx,9:13),1),1)];
                            elseif cs.nummed==12
                                L_Vi_Vphts = History{tk_Ci,2}(idx,3:8)+[History{tk_Ci,2}(idx,9:14),zeros(size(History{tk_Ci,2}(idx,9:14),1),1)];
                            end
                            W_phts = W0*exp(-sum(L_Vi_Vphts*muas',2));
                            
                            numerateur = zeros(1,6);
                            for ilayer=1:6
                                numerateur(1,ilayer) = sum(L_Vi_Vphts(:,ilayer).*W_phts,1);% sommation des photons
                            end
                            
                            PDP_Vi = numerateur/sum(W_phts);
                            PDPF_Vi(1:6,Ci) = PDP_Vi/chord;
                            
                            DPF = PDPF_Vi;
                        case 2
                            % Methode Michele
                            %: intégration à NIRS10 du code de Michèle
                            % Desjardins computeDirectProblem.m
                            if cs.nummed==6
                                photons_ppl = History{tk_Ci,2}(idx,3:size(History{tk_Ci,2},2));
                            elseif cs.nummed==11
                                photons_ppl = History{tk_Ci,2}(idx,3:8)+ History{tk_Ci,2}(idx,9:13);
                            end
                            % Il faut pondérer la moyenne sur les photons par le poids de chaque photon, donné par son atténuation dans le milieu. La probabilité qu'un photon ne soit pas absorbé dans un tissu est (exp(mua*parcours dans tissu)) (Hiraoka 1993).
                            photons_weight = exp(-1 * photons_ppl*muas');
                            
                            % Pour la paire pi (sm-dn), le DPF est la moyenne sur tous les
                            % photons, pondérée par leur atténuation, , de ce parcours
                            DPL(1:6,Ci) = photons_ppl * photons_weight ./ sum(photons_weight); % barycentre
                            DPF(1:6,Ci) = DPL(1,Ci)/(Cgp(1,Ci)*10);% ...divisée par la distance source-détecteur en mm                                     %%%%%%% Differential Pathlength Factor
                            
                            %                             % Parcours moyen dans la perturbation
                            %                             PPL(1,Ci) = photons_opl(6,Ci) * photons_weight(6,1) ./ sum(photons_weight(6,1));
                            %                             PPF(1,Ci) = PPL(1,Ci)/Cgp(1,Ci);
                            % Parcours moyen dans la zone BOLD
                            PPL(1,Ci) = Rphotons_ppl(7:11,Ci)'* photons_weight(1:5,1) ./ sum(photons_weight(1:5,1));
                            PPF(1,Ci) = PPL(1,Ci)/(Cgp(1,Ci)*10);
                            
                            PVF(1,Ci) = DPF(1,Ci)./PPF(1,Ci);% Partial pathlength factor and partial volume factor (PVF := DPF/PPF);
                    end
                    disp(['PVE calculated for channel : ' int2str(Ci)]);
                end
            catch % le fichier d histoire n a pas ete trouve...
                disp(['PVE failed for channel (must probable reason : source is not part of the ROI): ' int2str(Ci)]);
                PDPF_Vi(:,Ci) =-10;
            end
        end
        
        save(fullfile(cs_dir,'PDPF.mat'),'PDPF_Vi');
        NIRS.Dt.fir.PDPF = PDPF_Vi;
        %         save(fullfile(cs_dir,'DPF.mat'),'DPF');
        % NIRS.Dt.fir.DPF = DPF;
        %         save(fullfile(cs_dir,'PVE.mat'),'PVF');
        %         NIRS.Dt.fir.PVE = PVF;
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