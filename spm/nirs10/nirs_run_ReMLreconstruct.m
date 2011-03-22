function out = nirs_run_ReMLreconstruct(job)

load(job.NIRSmat{:});
% on recupere cs
dir_in = job.dir_in{:};
sep = strfind(dir_in,'\');
csn = dir_in(sep(end-1)+1:sep(end)-1);

itest=1;
while itest<length(NIRS.Cs.n)-1 && isempty(strfind(csn,NIRS.Cs.n{itest}))
    itest =itest+1;
end
i_cs =itest+1;
cs = NIRS.Cs.mcs{i_cs};

% volume 8bits et voxels isotropiques
% lecture du volume dans l'ordre de la matrice de sensitivite
fid=fopen(cs.b8i,'rb');
Yb8i = fread(fid);
fclose(fid);
Nvx = length(Yb8i);

% reperage des couches :
Yb8i_c1 = find(Yb8i==1);
Yb8i_c5 = find(Yb8i==5);

clear Yb8i;

% Le systeme qu'on resoud maintenant est le suivant :
% on est en un unique point temporel donc
%
%/ Y = X W betaWvx +epsilon_channel-noise
%\ betaWvx = omega
%
% bien piger si on est dqns le cas : 
% (la matrice de sensitivite est le FLUX....)
%/ Y = X W delta mua +epsilon_channel-noise
%\ delta mua = matrice de passage de mua a Hb. omega
%
% Y = un instant donne dans le fichier .nirs
% X = matrice sensitivite avec un point temporel
% W = matrice d'ondelettes spatialles
% betaWvx = concentrations dans les voxels
% epsilon_channel-noise = bruit dans les canaux
% omega = GOMBOUAMBA

%%%%%%% attention Huppert ne prend qu'une seule longueur d'onde !!!!
% Pairs....
Cmc = cs.C;
NC = length(Cmc); %Total number of measurements
Cwl=[];
Cwl = [Cwl NIRS.Cf.H.C.wl(Cmc)];
wl = unique(Cwl);

%% Y
fnirs = load(NIRS.Dt.fir.pp(1,4).p{:},'-mat');
t0 =3000;% on choisit au pif un point temporel
Y_t0 = fnirs.d(t0,Cmc)';

%% X
% p330
% on prend pour \Omegachapeau I en premiere approximation
Kmat = load(fullfile(dir_in,'sens.mat'));
Xdemi = Kmat.sens;
% maxi = max(Kmat.sens,[],2);
% for i=1:size(Xdemi,1)
%     Xdemi(i,:) = Xdemi(i,:)/maxi(i);
% end
% 1.2956 et 0.1940 sont les rapports entre les coefficients d'absortivite
% de HbO et HbR pour 690 et 830.
NCdemi = size(Xdemi,1)/2;
X = sparse([[Xdemi(1:NCdemi,:)*415.5 Xdemi(1:NCdemi,:)*2141.8];[Xdemi(NCdemi+1:end,:)*1008.0 Xdemi(NCdemi+1:end,:)*778.0]]);
% %Omega chapeau est le bruit d'observation. Ici on le prend egal a un
% Omegachapo = eye(NC);
% Kbarre = Omegachapo*K;
%% omega
beta_prior = zeros(2*Nvx,1);

%% W
% dans le cas d'une premiere reconstruction, on prend W =Id
W = sparse((1:2*Nvx),(1:2*Nvx),ones(1,2*Nvx)',2*Nvx,2*Nvx);

%% Qn : epsilon_channel-noise
for iwl=1:length(wl) %loop over number of wavelengths
    lst=find(Cwl==iwl); %List of all wavelength <idx>
    Qn{iwl}=sparse(lst,lst,ones(size(lst)),NC,NC);
end
%%  Qp : Covariance components for the parameters (4 total- 2 per HbO/HbR {layer 1; layer II})
% c1 : couche cortex : huppert 2
% c5 : couche skin   : huppert 1
%%% PREMIERE IDEE : on projette les ondelettes sur les couches en reperant
%%% quelle valeur appartient a quelle couche dans la longue matrice colonne
%%% (on fait comme si les voxels dans la matrice colonne etaient adjacents dans l'image)
%%% DEUXIEME IDEE
%- on calcule des wavelets dans un plan et on incline ce plan selon celui
%des sources et detecteurs

Qp{1}=sparse(Yb8i_c5,Yb8i_c5,zeros(length(Yb8i_c5),1),2*Nvx,2*Nvx);  %Skin layer- HbO
Qp{2}=sparse(Yb8i_c1,Yb8i_c1,ones(length(Yb8i_c1),1),2*Nvx,2*Nvx);  %Brain layer- HbO
Qp{3}=sparse(Yb8i_c5,Yb8i_c5,zeros(length(Yb8i_c5),1),2*Nvx,2*Nvx);  %Skin layer- HbR
Qp{4}=sparse(Yb8i_c1,Yb8i_c1,ones(length(Yb8i_c1),1),2*Nvx,2*Nvx);  %Brain layer- HbR

switch job.ReML_method
    case 0
        disp('code Huppert');
        meth = 'HUP';
        [lambda,beta_W,Stats]=nirs_run_DOT_REML(Y_t0,X*W',beta_prior,Qn,Qp);
        
        %Convert to the image domain and display
        beta = W'*beta_W;
        
    case 1
        disp('code spm_reml');
        meth = 'SPM';
        %Set up the extended covariance model by concatinating the measurement
        %and parameter noise terms
        Q=cell(length(Qn)+length(Qp),1);
    for idx=1:length(Qn)
        Q{idx}=blkdiag(Qn{idx},sparse(size(Qp{1},1),size(Qp{1},2))); % Build block diagonal matrix from Qn & Qp matrices
    end
    for idx2=1:length(Qp)
        Q{idx+idx2}=blkdiag(sparse(size(Qn{1},1),size(Qn{1},2)),Qp{idx2});
    end
        % sample covariance matrix Y*Y'
        YY = (Y_t0-mean(Y_t0))*(Y_t0-mean(Y_t0))';
        
        [C,h,Ph,F,Fa,Fc]=spm_reml(YY,X,Q);
        
        iC     = spm_inv(C);
        iCX    = iC*X;
        Cq = spm_inv(X'*iCX);
        
        beta = Cq*X'*iC*Y_t0;
    case 2
        disp('peudo inverse');
        meth = 'PInv';
        
        K=X;
        d=Y_t0;
        %from Edgar's code spm_lot_tikh
        %Tikhonov regularization parameters
        alpha = 1; %to be determined later
        f = (K*K' + alpha*eye(size(K,1))) \ d;
        beta = K'*f;
end
%%
% %Convert the Stats
% %%% juste les t stats
% Stats.tstat.Cov_beta = W'*Stats.tstat.Cov_beta*W;
% Stats.tstat.t=beta./sqrt(diag(Stats.tstat.Cov_beta));
% Stats.tstat.pval=2*tcdf(-abs(Stats.tstat.t),Stats.tstat.dfe);

V = spm_vol(cs.segR);
%Now, display the results
beta_4d = reshape(full(beta),[V.dim 2]);

beta_HbO = beta_4d(:,:,:,1);
beta_HbR = beta_4d(:,:,:,2);

% creation de nii :
V_O = struct('fname',fullfile(dir_in,['HbO_' meth '.nii']),...
    'dim',  V.dim,...
    'dt',   V.dt,...
    'pinfo',V.pinfo,...
    'mat',  V.mat);

V_O = spm_create_vol(V_O);
V_O = spm_write_vol(V_O, beta_HbO);

V_R = struct('fname',fullfile(dir_in,['HbR_' meth '.nii']),...
    'dim',  V.dim,...
    'dt',   V.dt,...
    'pinfo',V.pinfo,...
    'mat',  V.mat);

V_R = spm_create_vol(V_R);
V_R = spm_write_vol(V_R, beta_HbR);

% superpositions : on cree des images semi transparentes rouges ou bleues
% sur les anatomiques
spm_imcalc_ui({fullfile(dir_in,['HbO_' meth '.nii']);...
    cs.segR},...
    fullfile(dir_in,['HbO_anat_' meth '.nii']),...
    'i1+i2');

spm_imcalc_ui({fullfile(dir_in,['HbR_' meth '.nii']);...
    cs.segR},...
    fullfile(dir_in,['HbR_anat_' meth '.nii']),...
    'i1+i2');

out =1;
end