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
%_________________________________________________________________________
% Clément Bonnéry June 2011 : intégration à NIRS10 du code de Michèle
% Desjardins computeDirectProblem.m (pour tMCimg)

V0 = 0.05; % venous blood volume fraction
OEF0 = 1-0.6975; % baseline oxygen extraction fraction
HbR0 = 131.*(1-0.758); % baseline HbR concentration (=HbT0*(1-SO20)) in uM
TE = 30e-3; % 30 ms in our BOLD sequence
nu0 = 80.6; % s^-1, at 3T
alpha = -4.3 * nu0 * TE * V0 * OEF0 / HbR0; % Facteur de calibration du BOLD(relative change BOLD / uM) : BOLD(t)/BOLD(0) = alpha * [HbR(t)]

% Notations :
%  P = nombre de paires source-détecteur
%  V = nombre de voxels dans le volume de simulation Monte-Carlo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                                if tk_wl==690
                                    muas = [cs.par.gmPpties_l1(1,1),cs.par.wmPpties_l1(1,1),cs.par.csfPpties_l1(1,1),cs.par.skullPpties_l1(1,1),cs.par.scalpPpties_l1(1,1),cs.par.perturbationPpties_l1(1,1)];
                                elseif tk_wl ==830
                                    muas = [cs.par.gmPpties_l2(1,1),cs.par.wmPpties_l2(1,1),cs.par.csfPpties_l2(1,1),cs.par.skullPpties_l2(1,1),cs.par.scalpPpties_l2(1,1),cs.par.perturbationPpties_l2(1,1)];
                                end
                                [fluence, nPhoton] = lirehis_clm(tk_file,muas,1);
                                
                                
                            case 'mich'
                                %                 if sum(cs.par.perturbationPpties_l1 ==[0,0,0,0])/4 && sum(cs.par.perturbationPpties_l2 ==[0,0,0,0])/4
                                %                     ntissues = 5;
                                %                 else
                                %                     ntissues = 6;
                                %                 end
                                ntissues = 6;
                                [history] = lirehis_clm(t(k1,:),cs.par.nphotons,cs.NDkpt,ntissues,cs.par.numTimeGates);
                        end
                        History{tk_S} = history;
                    end
                end
                
                % cas caca, il faudra REFAIRE PLUS PROPRE
                Cid = NIRS.Cf.H.C.id;
                Cgp = NIRS.Cf.H.C.gp;
                for Ci=1:NIRS.Cf.H.C.N
                    try
                        %                        history_S_Ci = History{1,unique(Cid(2,Cid(1,:)==Ci))};
                        S_Ci = unique(Cid(2,Cid(1,:)== Ci));
                        D_Ci = unique(Cid(3,Cid(1,:)== Ci));
                        
                        idx = History{S_Ci}(:,1)==D_Ci;
                        photons_opl = sum(History{S_Ci}(idx,3:size(History{S_Ci},2)),2);% Parcours total (somme sur chaque tissu) de chaque photon issu de S et compté par D (chaque photon de la paire), en mm
                        %%%% pour d2tecteur le chemin dans le boldm
                        %%%% considerer le bold comme la perturbation
                        %%%% ///////
                        % % % %                         photons_opl_perturb = photons_opl(perturbation);
                        
                        % Il faut pondérer la moyenne sur les photons par le poids de chaque photon, donné par son atténuation dans le milieu. La probabilité qu'un photon ne soit pas absorbé dans un tissu est (exp(mua*parcours dans tissu)) (Hiraoka 1993).
                        photons_weight = exp(-1 * History{idxSrc}(idx,3:ntissus+2) * muaEachTiss');
                        % % % %                                 photons_weight_perturb = photons_weight(perturbation);
                        
                        % Pour la paire pi (sm-dn), le DPF est la moyenne sur tous les
                        % photons, pondérée par leur atténuation, , de ce parcours
                        if tk_wl==690
                            DPL(Ci,1) = photons_opl' * photons_weight ./ sum(photons_weight);
                            DPF(Ci,1) = DPL(Ci,1)./Cgp(Ci);% ...divisée par la distance source-détecteur en mm
                        elseif tk_wl ==830
                            DPL(Ci,2) = photons_opl' * photons_weight ./ sum(photons_weight);
                            DPF(Ci,2) = DPL(Ci,2)./Cgp(Ci);% ...divisée par la distance source-détecteur en mm
                        end
                        
                        % Parcours moyen dans la perturbation
                        PPL(pi,lambda) = parcoursPerturb' * weightPerturb ./ sum(weightPerturb);
                        PPF_lambda(pi) = PPL(pi,lambda)./(distSrcDet(sm,dn)*10);
                        
                        % Habituellement on considère un seul DPF pour une longueur d'onde et
                        % une distance source-détecteur données (moyenne sur toutes les paires,
                        % ça devrait se ressembler beaucoup si elles sont à la même distance)
                        DPFmean_lambda = mean(DPF_lambda);
                        DPFstd_lambda = std(DPF_lambda);
                        
                        % Partial pathlength factor and parital volume factor (PVF := DPF/PPF);
                        PPFmean_lambda = mean(PPF_lambda);
                        PPFstd_lambda = std(PPF_lambda);
                        PVF_lambda = DPF_lambda./PPF_lambda;
                        PVFmean_lambda = mean(PVF_lambda);
                        PVFstd_lambda = std(PVF_lambda);
                        
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