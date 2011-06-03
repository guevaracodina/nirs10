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
% Desjardins computeDirectProblem.m

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
% DelPreviousData  = job.DelPreviousData;
% try
%     NewNIRSdir = job.NewDirCopyNIRS.CreateNIRSCopy.NewNIRSdir;
%     NewDirCopyNIRS = 1;
% catch
%     NewDirCopyNIRS = 0;
% end
% Loop over subjects
for iSubj=1:size(job.NIRSmat,1)
    
    % Load NIRS.mat
    try
        NIRS = [];
        load(job.NIRSmat{iSubj,1});
        
        % gets current simulation cs
        sep = strfind(job.dir_in{:},'\');
        csn = dir_in(sep(end-1)+3:sep(end)-1);
        itest=1;
        while itest<length(NIRS.Cs.n) && (isempty(strfind(csn,NIRS.Cs.n{itest})) || length(csn)~=length(NIRS.Cs.n{itest}))
            itest =itest+1;
        end
        i_cs =itest;
        cs = NIRS.Cs.mcs{i_cs};
        
        switch cs.alg
            case 1 %MCX
            case 2 %tMCimg
                
                % ----------------------------------------------------------------------- %
                %% Géométrie : paires d'optodes
                
                
                % ICI : patch motrice ("moteur4x8") du CW6
                
                % Sources
                % nomsSrc = [{'s1-2'};
                %            {'s3-4'};
                %            {'s5-6'};
                %            {'s7-8'}];
                nomsSrc = [{'srcNo1'};
                    {'srcNo3'};
                    {'srcNo4'};
                    {'srcNo2'}];
                
                % Détecteurs correspondants - chaque ligne donne toutes les paires
                % associées à  une source
                % NE PAS CHANGER LA FORME "AAANOX" CAR CELA DÉTERMINE LE NO DE L'OPTODE
                % (NOM(6:END))
                nomsDet =  {[{'detNo9'}, {'detNo5'}, {'detNo8'}];
                    [{'detNo9'}, {'detNo12'}, {'detNo8'}, {'detNo6'}];
                    [{'detNo12'}, {'detNo11'}, {'detNo6'}, {'detNo7'}];
                    [{'detNo11'}, {'detNo10'}, {'detNo7'}]};
                % % Numéros de paires correspondants (informatif) PU BON:
                % pairnumbers = {[{'p1'}, {'p8'}, {'p9'}];
                %             [{'p2'}, {'p3'}, {'p10'}, {'p11'}];
                %             [{'p4'}, {'p5'}, {'p12'}, {'p13'}];
                %             [{'p6'}, {'p7'}, {'p14'}]};
                
                % Distance de chaque détecteur à sa source (dist. s-d pour chaque paire)
                % en cm (matrice paddée avec des 0)
                % (car unités des coeff. d'extinction (plus bas) sont cm^-1)
                % En utilisant les positions projetées à la surface du voume (réelles)
                distSrcDet = [ 3.2300    2.4800    2.2500    0; ...
                    3.2000    3.0700    2.6600    2.7500; ...
                    3.2000    3.1800    3.3100    2.7700; ...
                    4.0200    3.0600    2.8300    0;];
                
                
                %%
                %%% Calcul des phi et DPF pour chaque longueur d'onde... %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                lambdas = [690 830]; % CW6 : [690 830]. Sont inversées pour le CW5.
                
                % Paramètres
                nLambda = size(lambdas); % nb de longueurs d'onde
                P = 0;
                for src = 1:size(nomsDet,1)
                    P = P + size(nomsDet{src},2); % nb de paires
                end
                V = 0; % nb de voxels - sera lu, plus loin, dans les fichiers .cfg
                
                % Initialisations des matrices
                AA = [];
                % AA = zeros(2*P,2*V);         % Matrice de sensitivité
                % Sera initialisée plus loin, lorsqu'on a lu le nombre de voxels V dans les fichiers .cfg
                DPF = zeros(P,nLambda);        % Differential pathlength factor pour chaque paire et chaque longueur d'onde
                DPFmean= zeros(1,nLambda);     % Differential pathlength factor moyen pour chaque longueur d'onde
                DPL = zeros(P,nLambda);        % Parcours des photons pour chaque paire et chaque longueur d'onde
                % Note : DPF = DPL/distSrcDet;
                % Partial pathlength factors : for simulations with perturbation
                PPF = zeros(P,nLambda);        % Partial pathlength factor pour chaque paire et chaque longueur d'onde
                PPFmean= zeros(1,nLambda);     % Partial pathlength factor moyen pour chaque longueur d'onde
                PPL = zeros(P,nLambda);        % Parcours des photons dans la perturbation pour chaque paire et chaque longueur d'onde
                
                for lambda = 1:length(lambdas)
                    
                    %%% Liste des fichiers de sortie à lire %%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    % Demander à l'utilisateur : où sont les fichiers de sortie?
                    listFiles2pt = [];
                    listFileshis = [];
                    while isempty(listFiles2pt) || isempty(listFileshis)
                        dirMCoutputFiles = uigetdir('Z:\MonteCarloInVivo\P1S6\EntreesEtSortiesSimu\essai7_perturbe',...
                            'Dans quel dossier se trouvent les fichiers d''entrée de sortie de tMCimg?');
                        if (dirMCoutputFiles)
                            cd(dirMCoutputFiles);
                        else
                            break % If user cancels we get them out of the (otherwise infinite) while loop
                        end
                        
                        % Liste des fichiers .2pt et .his dans le dossier sélectionné
                        listFiles2pt = dir([dirMCoutputFiles '\*.2pt']);
                        listFileshis = dir([dirMCoutputFiles '\*.his']);
                        
                        %                 if ~isempty(listFiles2pt) && ~isempty(listFileshis)
                        %                     % listFiles est une structure nbFile x 1
                        %                     for i = 1:length(listFiles2pt)
                        %                         %nomsFichiers{i} = listFiles2pt(i).name(1:end-4);
                        %                         % On ne garde que le préfixe (enlève '_XXXnm.2pt')
                        %                         % Chaque nom de fichier est en paire, 2 par optode (pour 2
                        %                         % longueurs d'onde)
                        %                         if mod(i,2)==0 && mod(lambda,2)==0 % 830 nm
                        %                             nomsOptodes{ceil(i/2)} = listFiles2pt(i).name(1:end-10);
                        %                             nomsFichiers{ceil(i/2)} = listFiles2pt(i).name(1:end-4);
                        %                         elseif mod(i,2)==0 && mod(lambda,2)==1 % 690 nm
                        %                             nomsOptodes{ceil(i/2)} = listFiles2pt(i-1).name(1:end-10);
                        %                             nomsFichiers{ceil(i/2)} = listFiles2pt(i-1).name(1:end-4);
                        %                         end
                        %                     end
                        %
                        %                 else
                        %                     h=msgbox('Oups, pas de fichiers .2pt ou .his dans le dossier sélectionné!');
                        %                     uiwait(h);
                        %                     nomsFichiers = {};
                        %                 end
                        %             end
                        
                        
                        % ----------------------------------------------------------------------- %
                        
                        %%% Lire fichiers .his (durée et longueur du trajet de chaque photon)  %%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Lire 1 à 1 les fichiers .his
                        [allhistory nbPhotonsTotSimu ntissus muaEachTiss] = lirehis(nomsFichiers);
                        % allhis : cell contenant chaque history file
                        
                        % Pour chaque paire "pi" formée de la source "sm" et du détecteur "dn"...
                        pi = 0; % Numéro de la paire
                        phi0_lambda = zeros(P,1); % size : Px1
                        DPF_lambda = zeros(P,1); % size : Px1
                        PPF_lambda = zeros(P,1); % size : Px1
                        for sm = 1:size(nomsSrc,1)
                            for dn = 1:size(nomsDet{sm},2)
                                % paire pi : nomsSrc{sm} avec nomsDet{sm}{dn}
                                pi = pi+1;
                                % Index correspondant à sm dans allhis (déterminé par l'ordre
                                % dans lequel les noms de fichiers ont été passés par nomsFichiers)
                                idxSrc = find(strcmp(nomsSrc(sm),nomsOptodes)==1);
                                % % Index correspondant à dn dans allhis PAS PERTINENT
                                % idxDet = find(strcmp(nomsDet{sm}{dn},nomsOptodes)==1);
                                % Numéro d'optode du détecteur (pour la forme 'optodeNoX')
                                idxDet = str2num(nomsDet{sm}{dn}(6:end));
                                % idx : occurences de "détecteur dn ayant vu un photon
                                % de la source sm"
                                idx = find(allhistory{idxSrc}(:,1)==idxDet);
                                
                                %%% Calcul Diffenrential Pathlength Factor (DPF) %%%
                                %%% et Partial pathlength (PPF), le cas échéant  %%%
                                %%% (problème avec perturbation) %%%
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % Parcours total (somme sur chaque tissu) de chaque
                                % photon issu de sm et compté par dn (chaque photon de la
                                % paire), en mm
                                parcoursPhotons = sum(allhistory{idxSrc}(idx,3:ntissus+2),2);
                                % Parcours dans la perturbation (dernier tissu)
                                parcoursPerturb = allhistory{idxSrc}(idx,ntissus+2);
                                
                                % Il faut pondérer la moyenne sur les photons par le poids
                                % de chaque photon, donné par son atténuation dans le milieu.
                                % La probabilité qu'un photon ne soit pas absorbé dans un tissu
                                % est (exp(mua*parcours dans tissu)) (Hiraoka 1993).
                                weightPhotons = exp(-1 * allhistory{idxSrc}(idx,3:ntissus+2) * muaEachTiss');
                                weightPerturb = exp(-1 * allhistory{idxSrc}(idx,ntissus+2) * muaEachTiss(end));
                                
                                % Fluence au détecteur selon éq. (3) de Boas 2002
                                fluence(pi) = sum(weightPhotons)./nbPhotonsTotSimu;
                                
                                % Pour la paire pi (sm-dn), le DPF est la moyenne sur tous les
                                % photons, pondérée par leur atténuation, , de ce parcours
                                DPL(pi,lambda) = parcoursPhotons' * weightPhotons ./ sum(weightPhotons);
                                % ...divisée par la distance source-détecteur en mm
                                % (était en cm).
                                DPF_lambda(pi) = DPL(pi,lambda)./(distSrcDet(sm,dn)*10);
                                % Parcours moyen dans la perturbation
                                PPL(pi,lambda) = parcoursPerturb' * weightPerturb ./ sum(weightPerturb);
                                PPF_lambda(pi) = PPL(pi,lambda)./(distSrcDet(sm,dn)*10);
                                
                                disp('paire #:')
                                disp(pi)
                                disp(' # photons comptés sur cette paire :')
                                disp(length(idx))
                                disp(' ')
                                %disp(sum(weightPhotons))
                                
                                %%% Calcul phi0 %%% -> louche? Voir plus bas méthode plus
                                % fiable je crois...
                                % size : 1xP
                                % Ce calcul ne donne pas une valeur de phi0 normalisée comme
                                % le fait celui de phi à partir du fichier .2pt... par
                                % prudence je vais plutôt imiter le code de Boas
                                %%%%%%%%%%%%%%%%%%%
                                % phi0(pi) = phi0(sm,dn)
                                % = nb de photons issus de sm comptés par dn
                                nbPhotonsPair_lamda(pi) = length(idx);
                                %             phi0norm_lamda(pi) = phi0_lamda(pi) ./ nbPhotonsTotSimu;
                                %             % Comme on normalise phi dans lire2pt on veut aussi normaliser
                                %             % phi0
                                
                            end
                        end
                        
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
                        
                        % Libérer mémoire
                        %clear allhistory
                        
                        % ------------------------------------------------------------------- %
                        
                        %%% Lire fichiers .2pt (2-pt Green's functions ou densités de photons) %%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Lire 1 à 1 les fichiers .2pt (densités de photons), sommer sur le
                        % temps pour obtenir phi (on enregistre au passage cette matrice
                        % sommée qui prend moins d'espace mémoire)
                        [allPhi V volLims] = lire2pt(nomsFichiers);
                        [segVol volDimsxyz] = lirebin(nomsFichiers(1)); % ont tous le même fichier de segmentation donc on n'en lit qu'un
                        medium3D = segVol{1};
                        clear segVol
                        % allPhi : cell contenant chaque matrice phi (a été normalisé dans lire2pt)
                        % volLims : dimensions en voxels de la ROI dans laquelle on a calculé la simulation
                        % Format : [xmin ymin zmin; xmax ymax zmax]
                        % V : nombre de voxels dans la ROI de la simulation
                        % %%optodePositions : cell où chaque élément est un vecteur 1x3 de la
                        % %%position de la source dans le volume de simulation
                        % segVol : fichier de segmentation (volume segmenté en tissus)
                        
                        % Oups on ne veut pas réinitialiser à chaque longueur d'onde
                        %AA = zeros(2*P,2*V);         % Initialisation de la matrice de sensitivité
                        
                        %%% Calcul des éléments "L" (parcours effectif moyen) %%%
                        %%%               Matrice de sensitivité              %%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        % Pour chaque paire "pi" formée de la source "sm" et du détecteur "dn"...
                        pi = 0; % Numéro de la paire
                        % Matrice "L" (sensitivité)
                        L_lambda = zeros(P,V); % size : PxV
                        for sm = 1:size(nomsSrc,1)
                            for dn = 1:size(nomsDet{sm},2)
                                % paire pi : nomsSrc(sm) avec nomsDet(sm,dn)
                                pi = pi+1;
                                % Index correspondant à sm dans allphi
                                idxSrc = find(strcmp(nomsSrc(sm),nomsOptodes)==1);
                                % Index correspondant à dn dans allphi
                                idxDet = find(strcmp(nomsDet{sm}{dn},nomsOptodes)==1);
                                % Densités de photons phi
                                phi1 = allPhi{idxSrc}; % size nvx x nvy x nvz (V elements)
                                phi2 = allPhi{idxDet}; % idem
                                
                                %%% Calcul phi0 %%%
                                % size : 1xP
                                % Boas dans son code calcule son phi0 de cette façon, sans
                                % utiliser le fichier .his; il considère simplement la valeur
                                % de phi pour les photons issus de la source M dans le voxel
                                % correspondant à la position du détecteur N.
                                
                                % Lire les positions des optodes une fois projetées
                                fid = fopen('positionsProjected.txt');   % fichier que j'enregistre à la main
                                % Les positions doivent être dans l'ordre des optodes
                                % (1,2,3,...,12)
                                no_opt = 0;
                                while ~feof(fid)
                                    ligne = fgetl(fid);              % le lire ligne par ligne
                                    cel{1} = ligne;
                                    if isempty(ligne)
                                        break
                                    end
                                    no_opt = no_opt + 1;
                                    % Lire les positions des optodes
                                    positionyxz(no_opt,1:3) = floor(str2num(ligne(20:end-1)));
                                    % Vérifier qu'on est bien dans le tissu et non l'air
                                    if medium3D(positionyxz(no_opt,2),positionyxz(no_opt,1),positionyxz(no_opt,3))==0
                                        h = msgbox(['Attention la position trouvée du détecteur'...
                                            nomsDet{sm}{dn} 'est dans l''air!']);
                                        uiwait(h)
                                    end
                                end
                                fclose(fid);
                                
                                % Numéros d'optodes de sm et dn (dans la forme 'optodeNoX')
                                noDet = str2num(nomsDet{sm}{dn}(6:end));
                                noSrc = str2num(nomsSrc{sm}(6:end));
                                
                                % Phi0 pour la paire pi formée de la source sm et du détecteur
                                % dn est la valeur du phi de la source sm à la position du
                                % détecteur dn dans le volume. Ça pourrait aussi être la valeur
                                % du phi de la source dn à la position du détecteur sm. On
                                % prend la moyenne des deux.
                                phi0_lambdaS = phi1(positionyxz(noDet,2),...
                                    positionyxz(noDet,1),positionyxz(noDet,3));
                                phi0_lambdaD = phi2(positionyxz(noSrc,2),...
                                    positionyxz(noSrc,1),positionyxz(noSrc,3));
                                phi0_lambda(pi) = phi0_lambdaS/2 + phi0_lambdaD/2;
                                
                                phiProduct = phi1.*phi2; % idem (PxV)
                                L_lambda(pi,:) = phiProduct(:)' ./ phi0_lambda(pi); % size : 1xV
                                % Approximation de Rytov : le produit phi(s,v)*phi(v,d) est
                                % normalisé par phi0(s,d) - à l'opposé de l'approximation de
                                % Born (pas de normalisation)
                                
                                
                            end
                        end
                        
                        % Libérer mémoire
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
                        
                        % Pas nécessaire :
                        %phi0(:,lambda) = phi0_lambda;            % size : P x nLambda
                        
                    end
                    
                    % Enregistrer matrice de sensitivité pour l'avenir
                    save('MegaMatriceSensitivite.mat','AA');
                    clear L_lambda
                    save('allDPL','DPL','PPL','distSrcDet','nomsSrc','nomsDet','DPF','PVF');
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % ICI ON A GÉNÉRÉ LA MATRICE DE SENSITIVITÉ AA. LES LIGNES CI-DESSOUS FONT
                    % LA PROJECTION D'UN MASQUE BOLD À TRAVERS CETTE MATRICE DE SENSITIVITÉ,
                    % DANS L'ESPACE DES PAIRES SOURCES-DÉTECTEURS (PROBLÈME DIRECT).
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                    %%
                    %%% Matrice epsilons (coefficients d'extcintion) %%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Coefficients d'extinction selon Ted Huppert's HOMer (GetExtinctions.m)
                    % à, respectivement, 690 (ligne1, HbO HbR) et 830 (ligne2, HbO HbR)
                    
                    % [ Lamda1-HbO  Lambda1-HbR  (690=lambda1 avec le CW6)
                    %   Lambda2-HbO Lambda2-HbR] (830=lambda2 avec le CW6)
                    % size : nLambda x nHb (2x2)
                    
                    % source : http://omlc.ogi.edu/spectra/hemoglobin/summary.html (au 22
                    % juillet 2008)
                    if lambda(1)==830 && lambda(1)==690
                        coeff_ext = [ 974  693.04 ;   % POUR DONNÉES PRISES AVEC LE CW5
                            276  2051.96  ] .* 2.303 ./1e6; % en cm^-1 / (umol/L).
                    elseif lambda(1)==690 && lambda(1)==830
                        coeff_ext = [ 276  2051.96 ;   % POUR DONNÉES PRISES AVEC LE CW6
                            974  693.04  ] .* 2.303 ./1e6; % en cm^-1 / (umol/L).
                        % Le facteur 2.303 permet d'obtenir des coefficients d'absorption
                        % mu_a lorsqu'on mutliplie par la concentration en umol/L (toujours
                        % selon le site des données compilées par Scott Prahl!)
                    else
                        disp('Les coefficients utilisés dans ce script sont pour des longueurs d''onde de 690 et 830 nm!');
                    end
                    
                    
                    %% ----------------------------------------------------------------------- %
                    %%% Charger données IRMf %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Masque du changement de signal BOLD (rBOLD) en % par rapport au baseline,
                    % matrice en voxels, ici on ne garde que la ROI utilisée dans les
                    % simulations (volLims)
                    BOLDdata = load('maskBOLDP1SX.mat','-mat');
                    dPctfMRI_vox_BOLD = BOLDdata.BOLDmaskPct(volLims(1,1):volLims(2,1),...
                        volLims(1,2):volLims(2,2),volLims(1,3):volLims(2,3));
                    dPctfMRI_vox_BOLD(isnan(dPctfMRI_vox_BOLD)==1) = 0; % Il y a des NaN dans la matrice du masque
                    % car on a interpolé le scan BOLD aux dimensions, plus grandes, de l'IRM
                    % anatomique (même si ici on ne garde que les dimensions du MC). On les
                    % ramène à 0.
                    dConcfMRI_vox_HbR = dPctfMRI_vox_BOLD ./ alpha; % calibration du BOLD: BOLD=alpha*[HbR]
                    
                    % OUT OF MEMORY :
                    % dC = [dHbO; dHbR] = [0...0; dHbR];
                    %dConcfMRI_vox = zeros(2*V,1);
                    %dConcfMRI_vox(V+1:2*V) = dConcfMRI_vox_HbR(:);
                    %%% Projection des données IRMf en optique %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % % Référence : Huppert T., JBO 11(6), 2006
                    % %        dOD        =           A            *            E            *        dC
                    % %  (nLambda*P x 1)   (nLambda*P x nLambda*V)   (nLambda*V x nHbX*V)
                    % (nHbX*V x 1)
                    %dODproj_pairs = AA * kron(coeff_ext,eye(V)) * dConcfMRI_vox(:);
                    
                    % % Étape par étape :
                    % % dMua                  =        E            *        dC
                    % % (nLambda*V x 1)        (nLambda*V x nHbX*V)   (nHbX*V x 1)
                    % dMua_projHbR = kron(coeff_ext,eye(V)) * dConcfMRI_vox(:);
                    
                    % Ligne à ligne pour être memory-friendly:
                    
                    % Transformer dConc en dMua dans l'espace des voxels
                    % dMua      =        E            *        dC
                    % (1 x V)         (1 x nHbX)          (nHbX x V)
                    dConcfMRI_vox = zeros(2,V);
                    dConcfMRI_vox(2,:) = dConcfMRI_vox_HbR(:)'; % 1 x V
                    dMua_HbR_vox_lambda1 = coeff_ext(1,:) * dConcfMRI_vox; % 1 x V
                    dMua_HbR_vox_lambda2 = coeff_ext(2,:) * dConcfMRI_vox; % 1 x V
                    dMua_HbR_vox = [dMua_HbR_vox_lambda1'; dMua_HbR_vox_lambda2']; % (nHbX*V x 1)
                    
                    % Problème direct
                    % Projeter dMua des voxels dans l'espace des paires
                    % (nLambda*P x nLambda*V) * (nHbX*V x 1)
                    dODprojHbR_pairs = AA * dMua_HbR_vox;
                    % ----------------------------------------------------------------------- %
                    
                    
                    %% ----------------------------------------------------------------------- %
                    %%% Calcul des concentrations à partir des données optiques projetées %%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % % Loi de Beer-Lambert modifiée
                    
                    % % dOD(lambda,pi) = -ln(I/I0) = mua(lambda,pi)*dist. parcourue = mua*rSD*DPF
                    % % d'où :
                    % %     dmua =          dOD     ./     rSD    ./     DPF
                    % % (nLambda x P)  (nLambda x P) (nLambda x P)  (nLambda x P)
                    dmuaProj_pairs = reshape(dODprojHbR_pairs,[nLambda P]) ...
                        ./ repmat(distSrcDet(:)',[nLambda 1]) ...
                        ./ repmat(DPFmean(:),[1 P]) ;
                    % size : nLambda x P
                    % Ci-haut, on utilise un DPF moyen (pour toutes les paires) pour chaque
                    % longueur d'onde.
                    
                    % % Pour utiliser DPF propre à chaque paire:
                    % dmuaproj = reshape(dODproj,[nLambda p]) ...
                    %             ./ repmat(distSrcDet(:)',[nLambda 1]) ...
                    %             ./ DPF'; % (DPF était size PxnLambda)
                    
                    
                    % % Ensuite : dmua = coeff_ext * dConc
                    % % d'où :
                    % %     dConc     =    inv(coeff_ext)   *     dConc
                    % % (nLambda x P)   (nLambda x nLambda)   (nLambda x P)
                    dConcProj_pairs = inv(coeff_ext) * dmuaProj_pairs; % Il ne faut pas regarder HbO...
                    % HbR seulement (c'est la 2e ligne de l'éq. matricielle juste ci-haut)
                    dHbR_projBOLD = (coeff_ext(2,2) .* dmuaProj_pairs(1,:) - ...
                        coeff_ext(1,2) .* dmuaProj_pairs(2,:) ) ./ ...
                        (coeff_ext(2,2).*coeff_ext(1,1) - coeff_ext(2,1).*coeff_ext(1,2) );
                    
                    otherwise
                        disp('The algorithm with which the simulation has been runned is not recognised.')
                end
                
                % ----------------------------------------------------------------------- %
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
                    disp(['Conversion of optical intensities to hemoglobin ',...
                        'concentrations failed for subject ' int2str(Idx)]);
        end
    end
    out.NIRSmat = job.NIRSmat;
end