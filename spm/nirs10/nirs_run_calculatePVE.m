function out = nirs_run_calculatePVE(job)
% Calcul du facteur de volume partiel
%
% BOlDdata;  (SEULEMENT SI ON VEUT EFFECTUER UN PROBL�ME DIRECT, i.e. PROJETER CES DONN�ES DANS L'ESPACE DES PAIRES)
% Masque du changement de signal BOLD (rBOLD) en % par rapport au baseline,
% matrice en voxels (interpol�e � la r�solution de l'IRM anatomique)
% size : nHbX*V x 1
% il faut que dConcfMRI
% soit une matrice (d�roul�e en un vecteur) de V �l�ments (pour chaque
% esp�ce HbX concat�n�e) qui donne un masque qui vaut deltaHbO et R dans
% tout le volume.
% De l'analyse Neurolens : FIR estim�e normalis�e * valeur moyenne
% sur la ROI * masque qui vaut 1 dans la ROI et 0 ailleurs, somm�e (moyenn�e?) sur
% tmin:tmax.
%
% TRS : donne les donn�es n�cessaires au calcul du facteur de calibration du BOLD
% V0, OEF0 = (SO2in-SO2out)/SO2in = 1-ScO20, HbR0
%
% Quantit�s d'int�r�t calcul�es :
% 	DPF (differential pathlength factor);
% 	PPF (partial pathlength factor dans la perturbation, si on a roul� une simulation
% avec perturbation (on suppose que le dernier tissu est une perturbation)***);
% 	AA la matrice de sensitivit�
%
% *** Dans la simulation perturb�e, le volume (milieu de propagation) est comme suit :
% tissu 0=air, 1-5 = GM,WM,CSF,skull,scalp, tissu 6 = "perturbation" (avec propri�t�s
% optiques modifi�es � partir de celles de GM)
%_________________________________________________________________________
% Cl�ment Bonn�ry June 2011 : int�gration � NIRS10 du code de Mich�le
% Desjardins computeDirectProblem.m (pour tMCimg)

V0 = 0.05; % venous blood volume fraction
OEF0 = 1-0.6975; % baseline oxygen extraction fraction
HbR0 = 131.*(1-0.758); % baseline HbR concentration (=HbT0*(1-SO20)) in uM
TE = 30e-3; % 30 ms in our BOLD sequence
nu0 = 80.6; % s^-1, at 3T
alpha = -4.3 * nu0 * TE * V0 * OEF0 / HbR0; % Facteur de calibration du BOLD(relative change BOLD / uM) : BOLD(t)/BOLD(0) = alpha * [HbR(t)]

% Notations :
%  P = nombre de paires source-d�tecteur
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
                %%% Calcul des phi et DPF pour chaque longueur d'onde... %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                lambdas = [690 830]; % CW6 : [690 830]. Sont invers�es pour le CW5.
                               
                % Initialisations des matrices
%                 
%                 % Sera initialis�e plus loin, lorsqu'on a lu le nombre de voxels V dans les fichiers .cfg
%                 DPF = zeros(P,nLambda);        % Differential pathlength factor pour chaque paire et chaque longueur d'onde
%                 DPFmean= zeros(1,nLambda);     % Differential pathlength factor moyen pour chaque longueur d'onde
%                 DPL = zeros(P,nLambda);        % Parcours des photons pour chaque paire et chaque longueur d'onde
%                 % Note : DPF = DPL/distSrcDet;
%                 % Partial pathlength factors : for simulations with perturbation
%                 PPF = zeros(P,nLambda);        % Partial pathlength factor pour chaque paire et chaque longueur d'onde
%                 PPFmean= zeros(1,nLambda);     % Partial pathlength factor moyen pour chaque longueur d'onde
%                 PPL = zeros(P,nLambda);        % Parcours des photons dans la perturbation pour chaque paire et chaque longueur d'onde
                
                
                if sum(cs.par.perturbationPpties_l1 ==[0,0,0,0])/4 && sum(cs.par.perturbationPpties_l2 ==[0,0,0,0])/4
                    ntissues = 5;
                else
                    ntissues = 6;
                end
                
                for lambda = 1:length(lambdas)
                    
                    if strcmp(job.dir_in{:}(end),'\')
                        job.dir_in{:} = job.dir_in{:}(1:end-1);
                    end
                    [t,dummy] = spm_select('FPList',job.dir_in{:},'.his');
                    
                    for k1=1:size(t,1)
                        [allhistory nbPhotonsTotSimu ntissus muaEachTiss] = lirehis(t(k1,:),cs.par.nphotons,cs.NDkpt,ntissues,cs.par.numTimeGates);
                    end
                    
                    % Pour chaque paire "pi" form�e de la source "sm" et du d�tecteur "dn"...
                    pi = 0; % Num�ro de la paire
                    phi0_lambda = zeros(P,1); % size : Px1
                    DPF_lambda = zeros(P,1); % size : Px1
                    PPF_lambda = zeros(P,1); % size : Px1
                    for sm = 1:size(nomsSrc,1)
                        for dn = 1:size(nomsDet{sm},2)
                            % paire pi : nomsSrc{sm} avec nomsDet{sm}{dn}
                            pi = pi+1;
                            % Index correspondant � sm dans allhis (d�termin� par l'ordre
                            % dans lequel les noms de fichiers ont �t� pass�s par nomsFichiers)
                            idxSrc = find(strcmp(nomsSrc(sm),nomsOptodes)==1);
                            % % Index correspondant � dn dans allhis PAS PERTINENT
                            % idxDet = find(strcmp(nomsDet{sm}{dn},nomsOptodes)==1);
                            % Num�ro d'optode du d�tecteur (pour la forme 'optodeNoX')
                            idxDet = str2num(nomsDet{sm}{dn}(6:end));
                            % idx : occurences de "d�tecteur dn ayant vu un photon
                            % de la source sm"
                            idx = find(allhistory{idxSrc}(:,1)==idxDet);
                            
                            %%% Calcul Diffenrential Pathlength Factor (DPF) %%%
                            %%% et Partial pathlength (PPF), le cas �ch�ant  %%%
                            %%% (probl�me avec perturbation) %%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Parcours total (somme sur chaque tissu) de chaque
                            % photon issu de sm et compt� par dn (chaque photon de la
                            % paire), en mm
                            parcoursPhotons = sum(allhistory{idxSrc}(idx,3:ntissus+2),2);
                            % Parcours dans la perturbation (dernier tissu)
                            parcoursPerturb = allhistory{idxSrc}(idx,ntissus+2);
                            
                            % Il faut pond�rer la moyenne sur les photons par le poids
                            % de chaque photon, donn� par son att�nuation dans le milieu.
                            % La probabilit� qu'un photon ne soit pas absorb� dans un tissu
                            % est (exp(mua*parcours dans tissu)) (Hiraoka 1993).
                            weightPhotons = exp(-1 * allhistory{idxSrc}(idx,3:ntissus+2) * muaEachTiss');
                            weightPerturb = exp(-1 * allhistory{idxSrc}(idx,ntissus+2) * muaEachTiss(end));
                            
                            % Fluence au d�tecteur selon �q. (3) de Boas 2002
                            fluence(pi) = sum(weightPhotons)./nbPhotonsTotSimu;
                            
                            % Pour la paire pi (sm-dn), le DPF est la moyenne sur tous les
                            % photons, pond�r�e par leur att�nuation, , de ce parcours
                            DPL(pi,lambda) = parcoursPhotons' * weightPhotons ./ sum(weightPhotons);
                            % ...divis�e par la distance source-d�tecteur en mm
                            % (�tait en cm).
                            DPF_lambda(pi) = DPL(pi,lambda)./(distSrcDet(sm,dn)*10);
                            % Parcours moyen dans la perturbation
                            PPL(pi,lambda) = parcoursPerturb' * weightPerturb ./ sum(weightPerturb);
                            PPF_lambda(pi) = PPL(pi,lambda)./(distSrcDet(sm,dn)*10);
                            
                            disp('paire #:')
                            disp(pi)
                            disp(' # photons compt�s sur cette paire :')
                            disp(length(idx))
                            disp(' ')
                            %disp(sum(weightPhotons))
                            
                            %%% Calcul phi0 %%% -> louche? Voir plus bas m�thode plus
                            % fiable je crois...
                            % size : 1xP
                            % Ce calcul ne donne pas une valeur de phi0 normalis�e comme
                            % le fait celui de phi � partir du fichier .2pt... par
                            % prudence je vais plut�t imiter le code de Boas
                            %%%%%%%%%%%%%%%%%%%
                            % phi0(pi) = phi0(sm,dn)
                            % = nb de photons issus de sm compt�s par dn
                            nbPhotonsPair_lamda(pi) = length(idx);
                            %             phi0norm_lamda(pi) = phi0_lamda(pi) ./ nbPhotonsTotSimu;
                            %             % Comme on normalise phi dans lire2pt on veut aussi normaliser
                            %             % phi0
                            
                        end
                    end
                    
                    % Habituellement on consid�re un seul DPF pour une longueur d'onde et
                    % une distance source-d�tecteur donn�es (moyenne sur toutes les paires,
                    % �a devrait se ressembler beaucoup si elles sont � la m�me distance)
                    DPFmean_lambda = mean(DPF_lambda);
                    DPFstd_lambda = std(DPF_lambda);
                    
                    % Partial pathlength factor and parital volume factor (PVF := DPF/PPF);
                    PPFmean_lambda = mean(PPF_lambda);
                    PPFstd_lambda = std(PPF_lambda);
                    PVF_lambda = DPF_lambda./PPF_lambda;
                    PVFmean_lambda = mean(PVF_lambda);
                    PVFstd_lambda = std(PVF_lambda);
                    
                end     % Lib�rer m�moire
                %clear allhistory
                
                % Lib�rer m�moire
                clear allPhi phi1 phi2 phiProduct
                
                % ------------------------------------------------------------------- %
                
                % Matrices incluant toutes les longueurs d'onde %
                % --------------------------------------------- %
                DPF(:,lambda) = DPF_lambda;              % size : P x nLambda
                PVF(:,lambda) = PVF_lambda;              % size : P x nLambda
                DPFmean(lambda,1) = DPFmean_lambda;      % size : 1 x nLambda
                AA( (1+(P*(lambda-1))):(P+(P*(lambda-1))),...
                    (1+(V*(lambda-1))):(V+(V*(lambda-1))) ) = ...
                    L_lambda;       % size : 2P x 2V
                
                % Pas n�cessaire :
                %phi0(:,lambda) = phi0_lambda;            % size : P x nLambda
                
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                 end
                
                % Enregistrer matrice de sensitivit� pour l'avenir
                save('MegaMatriceSensitivite.mat','AA');
                clear L_lambda
                save('allDPL','DPL','PPL','distSrcDet','nomsSrc','nomsDet','DPF','PVF');
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % ICI ON A G�N�R� LA MATRICE DE SENSITIVIT� AA. LES LIGNES CI-DESSOUS FONT
                % LA PROJECTION D'UN MASQUE BOLD � TRAVERS CETTE MATRICE DE SENSITIVIT�,
                % DANS L'ESPACE DES PAIRES SOURCES-D�TECTEURS (PROBL�ME DIRECT).
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%
                %%% Matrice epsilons (coefficients d'extcintion) %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Coefficients d'extinction selon Ted Huppert's HOMer (GetExtinctions.m)
                % �, respectivement, 690 (ligne1, HbO HbR) et 830 (ligne2, HbO HbR)
                
                
                %% ----------------------------------------------------------------------- %
                %%% Charger donn�es IRMf %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Masque du changement de signal BOLD (rBOLD) en % par rapport au baseline,
                % matrice en voxels, ici on ne garde que la ROI utilis�e dans les
                % simulations (volLims)
                BOLDdata = load('maskBOLDP1SX.mat','-mat');
                dPctfMRI_vox_BOLD = BOLDdata.BOLDmaskPct(volLims(1,1):volLims(2,1),...
                    volLims(1,2):volLims(2,2),volLims(1,3):volLims(2,3));
                dPctfMRI_vox_BOLD(isnan(dPctfMRI_vox_BOLD)==1) = 0; % Il y a des NaN dans la matrice du masque
                % car on a interpol� le scan BOLD aux dimensions, plus grandes, de l'IRM
                % anatomique (m�me si ici on ne garde que les dimensions du MC). On les
                % ram�ne � 0.
                dConcfMRI_vox_HbR = dPctfMRI_vox_BOLD ./ alpha; % calibration du BOLD: BOLD=alpha*[HbR]
                
                % OUT OF MEMORY :
                % dC = [dHbO; dHbR] = [0...0; dHbR];
                %dConcfMRI_vox = zeros(2*V,1);
                %dConcfMRI_vox(V+1:2*V) = dConcfMRI_vox_HbR(:);
                %%% Projection des donn�es IRMf en optique %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % % R�f�rence : Huppert T., JBO 11(6), 2006
                % %        dOD        =           A            *            E            *        dC
                % %  (nLambda*P x 1)   (nLambda*P x nLambda*V)   (nLambda*V x nHbX*V)
                % (nHbX*V x 1)
                %dODproj_pairs = AA * kron(coeff_ext,eye(V)) * dConcfMRI_vox(:);
                
                % % �tape par �tape :
                % % dMua                  =        E            *        dC
                % % (nLambda*V x 1)        (nLambda*V x nHbX*V)   (nHbX*V x 1)
                % dMua_projHbR = kron(coeff_ext,eye(V)) * dConcfMRI_vox(:);
                
                % Ligne � ligne pour �tre memory-friendly:
                
                % Transformer dConc en dMua dans l'espace des voxels
                % dMua      =        E            *        dC
                % (1 x V)         (1 x nHbX)          (nHbX x V)
                dConcfMRI_vox = zeros(2,V);
                dConcfMRI_vox(2,:) = dConcfMRI_vox_HbR(:)'; % 1 x V
                dMua_HbR_vox_lambda1 = coeff_ext(1,:) * dConcfMRI_vox; % 1 x V
                dMua_HbR_vox_lambda2 = coeff_ext(2,:) * dConcfMRI_vox; % 1 x V
                dMua_HbR_vox = [dMua_HbR_vox_lambda1'; dMua_HbR_vox_lambda2']; % (nHbX*V x 1)
                
                % Probl�me direct
                % Projeter dMua des voxels dans l'espace des paires
                % (nLambda*P x nLambda*V) * (nHbX*V x 1)
                dODprojHbR_pairs = AA * dMua_HbR_vox;
                % ----------------------------------------------------------------------- %
                
                
                %% ----------------------------------------------------------------------- %
                %%% Calcul des concentrations � partir des donn�es optiques projet�es %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % % Loi de Beer-Lambert modifi�e
                
                % % dOD(lambda,pi) = -ln(I/I0) = mua(lambda,pi)*dist. parcourue = mua*rSD*DPF
                % % d'o� :
                % %     dmua =          dOD     ./     rSD    ./     DPF
                % % (nLambda x P)  (nLambda x P) (nLambda x P)  (nLambda x P)
                dmuaProj_pairs = reshape(dODprojHbR_pairs,[nLambda P]) ...
                    ./ repmat(distSrcDet(:)',[nLambda 1]) ...
                    ./ repmat(DPFmean(:),[1 P]) ;
                % size : nLambda x P
                % Ci-haut, on utilise un DPF moyen (pour toutes les paires) pour chaque
                % longueur d'onde.
                
                % % Pour utiliser DPF propre � chaque paire:
                % dmuaproj = reshape(dODproj,[nLambda p]) ...
                %             ./ repmat(distSrcDet(:)',[nLambda 1]) ...
                %             ./ DPF'; % (DPF �tait size PxnLambda)
                
                
                % % Ensuite : dmua = coeff_ext * dConc
                % % d'o� :
                % %     dConc     =    inv(coeff_ext)   *     dConc
                % % (nLambda x P)   (nLambda x nLambda)   (nLambda x P)
                dConcProj_pairs = inv(coeff_ext) * dmuaProj_pairs; % Il ne faut pas regarder HbO...
                % HbR seulement (c'est la 2e ligne de l'�q. matricielle juste ci-haut)
                dHbR_projBOLD = (coeff_ext(2,2) .* dmuaProj_pairs(1,:) - ...
                    coeff_ext(1,2) .* dmuaProj_pairs(2,:) ) ./ ...
                    (coeff_ext(2,2).*coeff_ext(1,1) - coeff_ext(2,1).*coeff_ext(1,2) );
                
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