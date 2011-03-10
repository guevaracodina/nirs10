function out = nirs_run_ReMLreconstruct(job)
% Demo of Restricted Maximum Likelihood estimation in DOT

load huppert_reml_demo.mat
%File contains:
%  SD-  Source-detector arrangement (see documentation for PMI toolbox from
%                   Harvard or HOMER software; www.nmr.mgh.harvard.edu/DOT)
%  X-   Optical forward model (including spectral priors) calculated from
%                   the PMI toolbox.  Y = X*[HbO; HbR]
%       NOTE- X is normalized by 1000 so both HbX and dOD are O(1-10)
%  Medium-  Structure describing the mesh used to generate the optical
%                   forward model (used here for display purposes)
%
%  SampleImage1-    Example of image (zero in layer 1)
%                   size = <x><y><z><{HbO,HbR}>
%                   size = 16 16 2  2
%
%  W -      the wavelet transform matrix 

%Run the code
clc;
disp('Tests on ReML code for the reconstruction of DOT images')
disp('Written by:  T. Huppert and Farras Abdelnour');
disp('Hijacked by:  C. Bonnery');
disp(' ');

%%Now, the actual data and reconstructions
load(job.NIRSmat{:});
% on recupere cs

% try ~sum(strcmp(NIRS.Cs.n,job.MCchoice.MC_nam))
%     csn = job.MC_nam;
%     else
%         csn = [job.MC_nam '1'];
% end
dir_in = job.dir_in{:};
sep = strfind(dir_in,'\');
csn = dir_in(sep(end-1)+1:sep(end)-1);

itest=1;
while itest<=length(NIRS.Cs.n) && isempty(strfind(csn,NIRS.Cs.n{itest}))
itest =itest+1;
end
i_cs =itest+1;

cs = NIRS.Cs.mcs{i_cs};

V = spm_vol(cs.segR);
Y = spm_read_vols(V);

%% Construct the covariance components and the hierarchical model
%Image prior 
% les concentrations sont sensees rester nulles
%%% est ce au4on cherche a etre sur un masque de la tete ???
% Beta_prior = zeros(size(donnees,2),1);  %Sparsity prior (akin to Tikhonov/Min Norm Est)
Beta_prior = zeros(size(Y));

% Construct the covariance components used in the model

%First, covariance components for the measurements
NC = NIRS.Cf.H.C.N;  %Total number of measurements
Cwl = NIRS.Cf.H.C.wl;
wl = unique(Cwl);
for iwl=1:length(wl) %loop over number of wavelengths
    lst=find(Cwl==iwl); %List of all wavelength <idx>
    Qn{iwl}=sparse(lst,lst,ones(size(lst)),NC,NC);
end

%Now, covariance components for the parameters (4 total- 2 per HbO/HbR {layer 1; layer II})

% c1 : couche cortex
% c5 : couche skin

Nvox=V.dim(1)*V.dim(2)*V.dim(3);
% on ne garde que deux couches dans lesquelles on reconstruit
lstskin=1:Nvox/2;  %This is a two-layer model
lstbrain=nvox/2+1:nvox;

%%
%This will act as a band-pass filter on each layer 
sigma_c5=6; 
sigma_c1=1;  %Sigma (see text) defines the attenuation at each frequency band
% If you increase sigma, this will act more as a low-pass filter (bias to
% low frequency).  If Skin> Brain, the bias will be to reconstruct the
% lower frequencies in layer 1


%This is a little messy since this demo does not provide all the code to
%calculate the original wavelets... so just 
% Number of stages
NS = 3;
temp = [ones(1,length(Medium.CompVol.X)/2^NS)];
for ii = NS:-1:1
    s1 = [temp ones(1,length(Medium.CompVol.X)/2^ii)/sigma_skin^(NS-ii+1)];
    temp = s1;
end
temp = [ones(1,length(Medium.CompVol.X)/2^NS)];
for ii = NS:-1:1
    s2 = [temp ones(1,length(Medium.CompVol.X)/2^ii)/sigma_brain^(NS-ii+1)];
    temp = s2;
end
%%% 1 skin
%%% 2 cortex
s1 = kron(s1',s1);  %I can do this as long as the image X/Y is square
s2 = kron(s2',s2);
skinWL_bias=s1(:);
brainWL_bias=s2(:);
%%

% Now, define the actual covariance components.  There are four; one for
% each layer times Hbo/HbR
%%%Note, we define this in the wavelet domain ...
Qp{1}=sparse(lstskin,lstskin,skinWL_bias,nvox*2,nvox*2);  %Skin layer- HbO
Qp{2}=sparse(lstbrain,lstbrain,brainWL_bias,nvox*2,nvox*2);  %Brain layer- HbO
Qp{3}=sparse(nvox+lstskin,nvox+lstskin,skinWL_bias,nvox*2,nvox*2);  %Skin layer- HbR
Qp{4}=sparse(nvox+lstbrain,nvox+lstbrain,brainWL_bias,nvox*2,nvox*2);  %Brain layer- HbR



%The actual model
disp('Computing multiple prior solution');
[lambda,Beta_W,Stats]=nirs_run_DOT_REML(Y,X*W',Beta_prior,Qn,Qp);
%lambda   - hyperparameters
%Beta_W   - the estimated image (in wavelet domain)
%Stats    - model Statistics (in the wavelet domain)

% %For comparison, the Tikhonov (equivelent) result
% disp('Computing Tikhonov/MNE solution');
% [lambda,Beta_W_Tikhonov,Stats]=nirs_run_DOT_REML(Y,X*W',Beta_prior,Qn,QpT);

%Convert to the image domain and display
%%% passage de beta a beta decompose sur l'espace des ondelettes de Daubechy
Beta = W'*Beta_W;
% Beta_Tikhonov = W'*Beta_W_Tikhonov;

%Convert the Stats
%%% juste les t stats
Stats.tstat.Cov_beta = W'*Stats.tstat.Cov_beta*W;
Stats.tstat.t=Beta./sqrt(diag(Stats.tstat.Cov_beta));
Stats.tstat.pval=2*tcdf(-abs(Stats.tstat.t),Stats.tstat.dfe);
%%%

%Now, display the results
%%% c'est biensur les beta qu'on affiche
Recon_Image=reshape(full(Beta),size(TrueImage));
% Recon_Image_Tikhonov=reshape(full(Beta_Tikhonov),size(TrueImage));

%%% Attention, TrueImage est  l'image en deux temps differents me
%%% semble-t-il...

%Scale all images the same
maxHbO1=max(max(abs([Recon_Image(:,:,1,1) TrueImage(:,:,1,1)])));
maxHbO2=max(max(abs([Recon_Image(:,:,2,1) TrueImage(:,:,2,1)])));
maxHbR1=max(max(abs([Recon_Image(:,:,1,2) TrueImage(:,:,1,2)])));
maxHbR2=max(max(abs([Recon_Image(:,:,2,2) TrueImage(:,:,2,2)])));

maxHbO1=max([maxHbO1 maxHbO2]); 
maxHbO2=maxHbO1;  %Remove this line if you don't want layer1 to be scaled the same as layer2
maxHbR1=max([maxHbR1 maxHbR2]); 
maxHbR2=maxHbR1;  %Remove this line if you don't want layer1 to be scaled the same as layer2


figure;
%%% on affiche les resultats qui ont ete calcules et ranges dans les
%%% matrices Recon_Image, Recon_Tikhonov_Image et TrueImage
%%% Medium.CompVol n'est que les echelles des distances selon X et Y
subplot(3,2,1); hold on;
imagesc(Medium.CompVol.X,Medium.CompVol.Y,squeeze(Recon_Image(:,:,1,1)));
for idx=1:size(SD.SrcPos,1); text(SD.SrcPos(idx,1),SD.SrcPos(idx,2),['S-' num2str(idx)]); end;
for idx=1:size(SD.DetPos,1); text(SD.DetPos(idx,1),SD.DetPos(idx,2),['D-' num2str(idx)]); end;
caxis([-maxHbO1 maxHbO1]); axis tight; axis off; colorbar;
title('REML Reconstructed Layer-I');

subplot(3,2,2); hold on;
imagesc(Medium.CompVol.X,Medium.CompVol.Y,squeeze(Recon_Image(:,:,2,1)));
for idx=1:size(SD.SrcPos,1); text(SD.SrcPos(idx,1),SD.SrcPos(idx,2),['S-' num2str(idx)]); end;
for idx=1:size(SD.DetPos,1); text(SD.DetPos(idx,1),SD.DetPos(idx,2),['D-' num2str(idx)]); end;
caxis([-maxHbO2 maxHbO2]); axis tight; axis off; colorbar;
title('REML Reconstructed Layer-II');

subplot(3,2,5); hold on;
imagesc(Medium.CompVol.X,Medium.CompVol.Y,squeeze(TrueImage(:,:,1,1)));
for idx=1:size(SD.SrcPos,1); text(SD.SrcPos(idx,1),SD.SrcPos(idx,2),['S-' num2str(idx)]); end;
for idx=1:size(SD.DetPos,1); text(SD.DetPos(idx,1),SD.DetPos(idx,2),['D-' num2str(idx)]); end;
caxis([-maxHbO1 maxHbO1]); axis tight; axis off; colorbar;
title('Truth Layer-I');

subplot(3,2,6); hold on;
imagesc(Medium.CompVol.X,Medium.CompVol.Y,squeeze(TrueImage(:,:,2,1)));
for idx=1:size(SD.SrcPos,1); text(SD.SrcPos(idx,1),SD.SrcPos(idx,2),['S-' num2str(idx)]); end;
for idx=1:size(SD.DetPos,1); text(SD.DetPos(idx,1),SD.DetPos(idx,2),['D-' num2str(idx)]); end;
caxis([-maxHbO2 maxHbO2]); axis tight; axis off; colorbar;
title('Truth Layer-II');
% end
end

% % % % % % % % % % % % % function out = nirs_run_ReMLreconstruct(job)
% % % % % % % % % % % % % % Achieve image segmentation after New Segment
% % % % % % % % % % % % % %_______________________________________________________________________
% % % % % % % % % % % % % % Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire
% % % % % % % % % % % % % % from Ted Huppert ReML code
% % % % % % % % % % % % % 
% % % % % % % % % % % % % % Clément Bonnéry
% % % % % % % % % % % % % % 2010-09
% % % % % % % % % % % % % 
% % % % % % % % % % % % % NIRS = job.NIRS;
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % NIRS.nirsfile = job.nirsfile
% % % % % % % % % % % % % % nirsfile = load('NIRS.nirsfile','-mat');
% % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % covariances inter and intra NIRS signals
% % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % intra (temporal covariance)
% % % % % % % % % % % % % % % % %-- sur OD : sur chaque paire
% % % % % % % % % % % % % % % % % inter (spatial covariance)
% % % % % % % % % % % % % % % % %-- on calcule les covariances entre les OD qui constituent les signaux
% % % % % % % % % % % % % % % % %directement enregistres sur les paires
% % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % %-- est ce qu'on devrait calculer la correlation entre HbO et HbR ????????
% % % % % % % % % % % % % % % % %(en prenant soin de correler HbO et -HbR disons...)
% % % % % % % % % % % % % 
% % % % % % % % % % % % % n_pairs = size(NIRS.nirs_file.d,2)/2;
% % % % % % % % % % % % % % SAME wavelength
% % % % % % % % % % % % % for n_wl = 1:2
% % % % % % % % % % % % %     for i = 1:n_pairs
% % % % % % % % % % % % %         for j = 1:n_pairs
% % % % % % % % % % % % %             Q{(n_wl-1)*n_pairs+i,(n_wl-1)*n_pairs+j} = xcorr(NIRS.nirs_file.d(:,(n_wl-1)*n_pairs+i),NIRS.nirs_file.d(:,(n_wl-1)*n_pairs+j));
% % % % % % % % % % % % %             Q{i,(n_wl-1)*n_pairs+j} = xcorr(NIRS.nirs_file.d(:,i),NIRS.nirs_file.d(:,n_pairs+j));
% % % % % % % % % % % % %             Q{(n_wl-1)*n_pairs+i,j} = xcorr(NIRS.nirs_file.d(:,n_pairs+i),NIRS.nirs_file.d(:,j));
% % % % % % % % % % % % %         end
% % % % % % % % % % % % %     end
% % % % % % % % % % % % % end
% % % % % % % % % % % % % 
% % % % % % % % % % % % % [C,h,Ph,F,Fa,Fc] = spm_reml_sc(YY,X,Q,NIRS.nirs_file.size(d,1),-32,256,V);
% % % % % % % % % % % % % % ReML estimation of covariance components from y*y' - proper components
% % % % % % % % % % % % % % FORMAT [C,h,Ph,F,Fa,Fc] = spm_reml_sc(YY,X,Q,N,[hE,hC,V]);
% % % % % % % % % % % % % %
% % % % % % % % % % % % % % YY  - (m x m) sample covariance matrix Y*Y'  {Y = (m x N) data matrix}
% % % % % % % % % % % % % % X   - (m x p) design matrix
% % % % % % % % % % % % % % Q   - {1 x q} covariance components
% % % % % % % % % % % % % % N   - number of samples
% % % % % % % % % % % % % %
% % % % % % % % % % % % % % hE  - hyperprior expectation in log-space [default = -32]
% % % % % % % % % % % % % % hC  - hyperprior covariance  in log-space [default = 256]
% % % % % % % % % % % % % % V   - fixed covariance component
% % % % % % % % % % % % % %
% % % % % % % % % % % % % % C   - (m x m) estimated errors = h(1)*Q{1} + h(2)*Q{2} + ...
% % % % % % % % % % % % % % h   - (q x 1) ReML hyperparameters h
% % % % % % % % % % % % % % Ph  - (q x q) conditional precision of log(h)
% % % % % % % % % % % % % %
% % % % % % % % % % % % % % F   - [-ve] free energy F = log evidence = p(Y|X,Q) = ReML objective
% % % % % % % % % % % % % %
% % % % % % % % % % % % % % Fa  - accuracy
% % % % % % % % % % % % % % Fc  - complexity (F = Fa - Fc)
% % % % % % % % % % % % % %
% % % % % % % % % % % % % % Performs a Fisher-Scoring ascent on F to find MAP variance parameter
% % % % % % % % % % % % % % estimates.  NB: uses weakly informative log-normal hyperpriors.
% % % % % % % % % % % % % % See also spm_reml for an unconstrained version that allows for negative
% % % % % % % % % % % % % % hyperparameters
% % % % % % % % % % % % % %
% % % % % % % % % % % % % %__________________________________________________________________________
% % % % % % % % % % % % % %      spm_reml_sc: positivity constraints on covariance parameters
% % % % % % % % % % % % % 
% % % % % % % % % % % % % % on resout à t donné
% % % % % % % % % % % % % 
% % % % % % % % % % % % % %get anatomical and functional datas
% % % % % % % % % % % % % %-> anatomical datas : 5 layer segmented image
% % % % % % % % % % % % % 
% % % % % % % % % % % % % 
% % % % % % % % % % % % % 
% % % % % % % % % % % % % %-> functional datas :
% % % % % % % % % % % % % %   -- position of sources and detectors
% % % % % % % % % % % % % %   -- HbO and HbT [in SD pairs domain]
% % % % % % % % % % % % % %   -- HbO and HbT hemodynamic response ?? (prior temporel...)
% % % % % % % % % % % % % %   -- ??
% % % % % % % % % % % % % out{1} =1;
% % % % % % % % % % % % % return
% % % % % % % % % % % % % 
% % % % % % % % % % % % % % function out = nirs_run_reconstruction(job)
% % % % % % % % % % % % % % % Achieve image segmentation after New Segment
% % % % % % % % % % % % % % %_______________________________________________________________________
% % % % % % % % % % % % % % % Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire
% % % % % % % % % % % % % % % from Ted Huppert ReML code
% % % % % % % % % % % % % %
% % % % % % % % % % % % % % % Clément Bonnéry
% % % % % % % % % % % % % % % 2010-09
% % % % % % % % % % % % % %
% % % % % % % % % % % % % %
% % % % % % % % % % % % % % if(~exist('K'))
% % % % % % % % % % % % % %     K=35;  %Max # of iterations of REML code
% % % % % % % % % % % % % % end
% % % % % % % % % % % % % %
% % % % % % % % % % % % % % % on resout à t donné
% % % % % % % % % % % % % %
% % % % % % % % % % % % % % %get anatomical and functional datas
% % % % % % % % % % % % % % %-> anatomical datas : 5 layer segmented image
% % % % % % % % % % % % % % %-> functional datas :
% % % % % % % % % % % % % % %   -- position of sources and detectors
% % % % % % % % % % % % % % %   -- HbO and HbT [in SD pairs domain]
% % % % % % % % % % % % % % %   -- HbO and HbT hemodynamic response ?? (prior temporel...)
% % % % % % % % % % % % % % %   -- ??
% % % % % % % % % % % % % %
% % % % % % % % % % % % % % %set the model (write down explicitly the system)
% % % % % % % % % % % % % % % From SPM book p 283______________________________________________________
% % % % % % % % % % % % % %
% % % % % % % % % % % % % % % setting of the model : augmentation to embody priors in error covariance
% % % % % % % % % % % % % % X = [X,X;speye(size(X,1))];
% % % % % % % % % % % % % % y = [y;zeros(size(X,2));eta_beta];
% % % % % % % % % % % % % % C_beta = sparse(i,j,C_beta,m,n);
% % % % % % % % % % % % % % C_epsilon = C_beta;
% % % % % % % % % % % % % %
% % % % % % % % % % % % % % lambda = ones(length(Q),1);  %Initial guess of lambda.
% % % % % % % % % % % % % % tol = 1E-4; %REML goes till tolerance (or max iter)
% % % % % % % % % % % % % % t = 256;
% % % % % % % % % % % % % % dF = Inf;   %Initial decent
% % % % % % % % % % % % % % cnt = 0;    %This is a bookkeeping param for display purposes
% % % % % % % % % % % % % %
% % % % % % % % % % % % % %
% % % % % % % % % % % % % % % building Q matrixes______________________________________________________
% % % % % % % % % % % % % %
% % % % % % % % % % % % % %
% % % % % % % % % % % % % % % until convergence EM algorithm___________________________________________
% % % % % % % % % % % % % % for iter=1:K
% % % % % % % % % % % % % %
% % % % % % % % % % % % % %     % E-Step :
% % % % % % % % % % % % % %     %______________________________________________________________________
% % % % % % % % % % % % % %     %   -- C_epsilon
% % % % % % % % % % % % % %     for i = 1:length(Q)
% % % % % % % % % % % % % %         C_epsilon = C_epsilon + Q{i}*lambda(i);
% % % % % % % % % % % % % %     end
% % % % % % % % % % % % % %     %   -- C_beta|y
% % % % % % % % % % % % % %     iC_epsilon = inv(C_epsilon);
% % % % % % % % % % % % % %     Xt_iC_epsilon = X' * iC_epsilon;
% % % % % % % % % % % % % %     Xt_iC_epsilon_X = Xt_iC_epsilon * X;
% % % % % % % % % % % % % %     C_beta_y = inv(Xt_iC_epsilon_X);  %Estimate of covariance of beta given the measurements
% % % % % % % % % % % % % %     % -- n_beta_y : jamais utilise...........
% % % % % % % % % % % % % %
% % % % % % % % % % % % % %     % M-Step :
% % % % % % % % % % % % % %     %______________________________________________________________________
% % % % % % % % % % % % % %     %   -- P
% % % % % % % % % % % % % %     P = iC_epsilon - (iC_epsilon*X)*C_beta_y*Xt_iC_epsilon;
% % % % % % % % % % % % % %     %   -- g_i and H_ij
% % % % % % % % % % % % % %     PY=P*Y;
% % % % % % % % % % % % % %     for i=1:size(Q,1)
% % % % % % % % % % % % % %         PQ_i{i}=P*Q{i};
% % % % % % % % % % % % % %     end
% % % % % % % % % % % % % %
% % % % % % % % % % % % % %     for i = 1:size(Q,1)
% % % % % % % % % % % % % %         PQ = PQ_i{i};
% % % % % % % % % % % % % %         PQt=PQ';
% % % % % % % % % % % % % %         g(i,1) = -0.5*trace(PQ)*exp(lambda(i)) + 0.5*PY'*Q{i}*PY*exp(lambda(i));
% % % % % % % % % % % % % %         for j = i:size(Q,1)
% % % % % % % % % % % % % %             PQj = PQ_i{j};
% % % % % % % % % % % % % %             H(i,j) = -0.5*sum(sum(PQt.*PQj))*exp(lambda(i)+lambda(j));
% % % % % % % % % % % % % %             H(j,i)=H(i,j);
% % % % % % % % % % % % % %         end
% % % % % % % % % % % % % %     end
% % % % % % % % % % % % % %
% % % % % % % % % % % % % %     %Now update the lambda.  dLambda = -inv(H)*g
% % % % % % % % % % % % % %     %     I=eye(size(H,1));
% % % % % % % % % % % % % %     dlambda = -inv(H)*g;%(expm(H*t) - I)*inv(H)*g;
% % % % % % % % % % % % % %     lambda = lambda + dlambda;
% % % % % % % % % % % % % %
% % % % % % % % % % % % % %     df    = g'*dlambda;
% % % % % % % % % % % % % %     if df > dF - exp(-4), t = max(2,t/2); end %retune the regularization if req.????????????????
% % % % % % % % % % % % % %     dF    = df;
% % % % % % % % % % % % % %
% % % % % % % % % % % % % %     for c=1:cnt, fprintf('\b'); end
% % % % % % % % % % % % % %     cnt=fprintf('%-5s: %i %5s%e','  ReML Iteration',iter,'...',full(dF));
% % % % % % % % % % % % % %
% % % % % % % % % % % % % %     if dF < tol, break; end
% % % % % % % % % % % % % %
% % % % % % % % % % % % % % end
% % % % % % % % % % % % % %
% % % % % % % % % % % % % % %Now, put the final pieces together
% % % % % % % % % % % % % % Beta= C_beta_y * Xt_iCe * Y;
% % % % % % % % % % % % % %
% % % % % % % % % % % % % % fprintf('\n');
% % % % % % % % % % % % % % out{1} =1;
% % % % % % % % % % % % % % return