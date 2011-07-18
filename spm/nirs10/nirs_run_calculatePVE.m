function out = nirs_run_calculatePVE(job)
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
% 	AA la matrice de sensitivité
%
% *** Dans la simulation perturbée, le volume (milieu de propagation) est comme suit :
% tissu 0=air, 1-5 = GM,WM,CSF,skull,scalp, tissu 6 = "perturbation" (avec propriétés
% optiques modifiées à partir de celles de GM)
%
% V0 = 0.05; % venous blood volume fraction
% OEF0 = 1-0.6975; % baseline oxygen extraction fraction
% HbR0 = 131.*(1-0.758); % baseline HbR concentration (=HbT0*(1-SO20)) in uM
% TE = 30e-3; % 30 ms in our BOLD sequence
% nu0 = 80.6; % s^-1, at 3T
% alpha = -4.3 * nu0 * TE * V0 * OEF0 / HbR0; % Facteur de calibration du BOLD(relative change BOLD / uM) : BOLD(t)/BOLD(0) = alpha * [HbR(t)]
%
%_________________________________________________________________________
% Clément Bonnéry June 2011 : intégration à NIRS10 du code de Michèle
% Desjardins computeDirectProblem.m (pour tMCimg)

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
        
        if isempty(job.dir_in{:})
            [dir_in, dummy1, dummy2] = fileparts(job.NIRSmat{iSubj,1});
        else
            dir_in = job.dir_in{:};
        end
        
        % gets current simulation cs
        sep = strfind(dir_in,'\');
        csn = dir_in(sep(end-1)+3:sep(end)-1);
        itest=1;
        while itest<length(NIRS.Cs.n) && (isempty(strfind(csn,NIRS.Cs.n{itest})) || length(csn)~=length(NIRS.Cs.n{itest}))
            itest =itest+1;
        end
        i_cs =itest;
        cs = NIRS.Cs.mcs{i_cs};
        
        switch cs.alg
            case 1 %MCX
                [t,dummy] = spm_select('FPList',fullfile(cs.p,cs.dir),'.mch');
                for i=1:size(t,1)
                    [data, header]=loadmch(t(i,:));
                end
                
                
                
            case 2 %tMCimg
                if strcmp(job.dir_in{:}(end),'\')
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
                        
                        if tk_wl==690
                            muas = [cs.par.gmPpties_l1(1,1),cs.par.wmPpties_l1(1,1),cs.par.csfPpties_l1(1,1),cs.par.skullPpties_l1(1,1),cs.par.scalpPpties_l1(1,1),cs.par.perturbationPpties_l1(1,1)];
                        elseif tk_wl ==830
                            muas = [cs.par.gmPpties_l2(1,1),cs.par.wmPpties_l2(1,1),cs.par.csfPpties_l2(1,1),cs.par.skullPpties_l2(1,1),cs.par.scalpPpties_l2(1,1),cs.par.perturbationPpties_l2(1,1)];
                        end
                        
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
                                [history] = lirehis_clm(t(k1,:),cs.par.nphotons,cs.NDkpt,ntissues,cs.par.numTimeGates);
                                count = count+1;
                        end
                        History{count,1} = tk_file;
                        History{count,2} = history;
                    end
                end
                
                % cas caca, il faudra REFAIRE PLUS PROPRE
                Cid = NIRS.Cf.H.C.id;
                Cgp = NIRS.Cf.H.C.gp;
                for Ci=1:NIRS.Cf.H.C.N
                    try
                        S_Ci = unique(Cid(2,Cid(1,:)== Ci));
                        D_Ci = unique(Cid(3,Cid(1,:)== Ci));
                        
                        for iwl =1:length(NIRS.Cf.dev.wl)
                            %%% on reconstruit le nom
                            if S_Ci < 10, S_Cin = ['0' int2str(S_Ci)]; else S_Cin = int2str(S_Ci);end
                            tk_n = ['S_No' S_Cin '_' int2str(NIRS.Cf.dev.wl(iwl)) 'nm'];
                            i=1;
                            while strcmp(tk_n,History{i,1})==0
                                i=i+1;
                            end
                            tk_Ci = i;% indice de la bonne ligne pour le canal Ci
                            
                            % maintenant les detecteurs sont ranges comme
                            % dans le fichier de config, il y en a le
                            % nombre de P-1 et c est dans lordre en
                            % enlevant la bonne source
                            %%% on veut D_Ci
                            D_Ci_cfg = 7+D_Ci;
                            
                            idx = History{tk_Ci,2}(:,1)==D_Ci_cfg;
                            if sum(idx)==0
                                if tk_wl==690
                                    PVF(1,Ci)=0;
                                else
                                    PVF(2,Ci)=0;
                                end
                            else
                                photons_opl(:,Ci) = sum(History{tk_Ci,2}(idx,3:size(History{tk_Ci,2},2)),1);% Parcours total (somme sur chaque tissu) de chaque photon issu de S et compté par D (chaque photon de la paire), en mm
                                
                                % Il faut pondérer la moyenne sur les photons par le poids de chaque photon, donné par son atténuation dans le milieu. La probabilité qu'un photon ne soit pas absorbé dans un tissu est (exp(mua*parcours dans tissu)) (Hiraoka 1993).
                                photons_weight = exp(-1 * photons_opl(:,Ci).*muas');
                                
                                % Pour la paire pi (sm-dn), le DPF est la moyenne sur tous les
                                % photons, pondérée par leur atténuation, , de ce parcours
                                if tk_wl==690
                                    DPL(1,Ci) = photons_opl'* photons_weight ./ sum(photons_weight);
                                    DPF(1,Ci) = DPL(1,Ci)/Cgp(1,Ci);% ...divisée par la distance source-détecteur en mm
                                    
                                    % Parcours moyen dans la perturbation
                                    PPL(1,Ci) = photons_opl(6,1) * photons_weight(6,1) ./ sum(photons_weight(6,1));
                                    PPF(1,Ci) = PPL(1,Ci)/Cgp(1,Ci);
                                    
                                    PVF(1,Ci) = DPF(1,Ci)./PPF(1,Ci);% Partial pathlength factor and partial volume factor (PVF := DPF/PPF);
                                    
                                elseif tk_wl ==830
                                    DPL(2,Ci) = photons_opl'* photons_weight ./ sum(photons_weight);
                                    DPF(2,Ci) = DPL(2,Ci)/Cgp(1,Ci);% ...divisée par la distance source-détecteur en mm
                                    
                                    % Parcours moyen dans la perturbation
                                    PPL(2,Ci) = photons_opl(6,1) * photons_weight(6,1) ./ sum(photons_weight(6,1));
                                    PPF(2,Ci) = PPL(2,Ci)/Cgp(1,Ci);
                                    
                                    PVF(2,Ci) = DPF(2,Ci)./PPF(2,Ci);% Partial pathlength factor and partial volume factor (PVF := DPF/PPF);
                                end
                            end
                        end
                        
                    catch % le fichier d histoire n a pas ete trouve...
                        disp(['PVE failed for channel : ' Ci]);
                    end
                end
            otherwise
                disp('The algorithm with which the simulation has been runned is not recognised.')
        end
        
        % ----------------------------------------------------------------------- %
        %[dir1,fil1,ext1] = fileparts(NIRS.Dt.s.p);
        if NewDirCopyNIRS
            dir2 = [dir1 filesep NewNIRSdir];
            if ~exist(dir2,'dir'), mkdir(dir2); end;
            outfile = fullfile(dir2,[prefix fil1 ext1]);
        else
            outfile = fullfile(dir1,[prefix fil1 ext1]);
        end
        if DelPreviousData
            delete(rDtp{f,1});
        end
        if NewDirCopyNIRS
            newNIRSlocation = fullfile(dir2,'NIRS.mat');
            save(newNIRSlocation,'NIRS');
            job.NIRSmat{Idx,1} = newNIRSlocation;
        else
            save(job.NIRSmat{Idx,1},'NIRS');
        end
    catch
        disp(['Calculus of PVE failed for subject ' int2str(iSubj)]);
    end
end
out.NIRSmat = job.NIRSmat;
end