% Demo of Restricted Maximum Likelihood estimation in DOT
% as described in:
%
% (in review to ???)

% Ted Huppert
% integrqted in nirs10 by clement

load huppert_reml_demo.mat
%File contains:
%  SD-  Source-detector arrangement (see documentation for PMI toolbox from
%                   Harvard or HOMER software; www.nmr.mgh.harvard.edu/DOT)
%  X-   Optical forward model (including spectral priors) calculated from
%                   the PMI toolbox.  Y = X*[HbO; HbR]
%       NOTE- X is normalized by 1000 so both HbX and dOD are O(1-10)
%  Medium-  Structure describing the mesh used to generate the optical
%                   forward model (used here for display purposes)
%  SampleImage1-    Example of image (zero in layer 1)
%                   size = <x><y><z><{HbO,HbR}>
%
%  SampleImage2-    Example of image (non-zero in layer 1)
%
%  W -      the wavelet transform matrix 

%Note- this example is a bit smaller than the one described in the paper
%because I want this to run faster.  The model used in the paper took
%~20min per fit.


CNR=5;  %Define the contrast-to-noise of the data


%Run the code
clc;
disp('Example of ReML code for the reconstruction of DOT images')
disp('Written by:  T. Huppert and Farras Abdelnour');
disp('University of Pittsburgh');
disp(' ');

%% 
%Construct the covariance components and the heirachical model

%Image prior 
Beta_prior = zeros(size(X,2),1);  %Sparsity prior (akin to Tikhonov/Min Norm Est)


%% Construct the covariance components used in the model

%First, covariance components for the measurements
ny=size(SD.MeasList,1);  %Total number of measurements
for idx=1:length(SD.Lambda) %loop over number of wavelengths
    lst=find(SD.MeasList(:,4)==idx); %List of all wavelength <idx>
    Qn{idx}=sparse(lst,lst,ones(size(lst)),ny,ny);
end


%Now, covariance components for the parameters (4 total- 2 per HbO/HbR {layer 1; layer II})
nvox=Medium.nVox;
lstskin=1:nvox/2;  %This is a two-layer model
lstbrain=nvox/2+1:nvox;

%%
%This will act as a band-pass filter on each layer 
sigma_skin=6; 
sigma_brain=1;  %Sigma (see text) defines the attenuation at each frequency band
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
s1 = kron(s1',s1);  %I can do this as long as the image X/Y is square
s2 = kron(s2',s2);
skinWL_bias=s1(:);
brainWL_bias=s2(:);
%%

% Now, define the actual covariance components.  There are four; one for
% each layer times Hbo/HbR
%Note, we define this in the wavelet domain ...
Qp{1}=sparse(lstskin,lstskin,skinWL_bias,nvox*2,nvox*2);  %Skin layer- HbO
Qp{2}=sparse(lstbrain,lstbrain,brainWL_bias,nvox*2,nvox*2);  %Brain layer- HbO
Qp{3}=sparse(nvox+lstskin,nvox+lstskin,skinWL_bias,nvox*2,nvox*2);  %Skin layer- HbR
Qp{4}=sparse(nvox+lstbrain,nvox+lstbrain,brainWL_bias,nvox*2,nvox*2);  %Brain layer- HbR


% This is a single prior (MNE) which will give us results equvelent to the
% MNE prior Beta = inv(X'*X+lambda*I)*X'*Y
QpT{1}=speye(nvox*2);  %For equivelent to Tikhonov/MNE;


%%Now, the actual data and reconstructions

for idx=1:2
    %Two passes.  
    Contrast=std(X * SampleImage1(:));
    
    if(idx==1)
    %Example 1: No noise in skin layer
        disp('Example I- no superficial noise');
        TrueImage=SampleImage1;
    elseif(idx==2)
        disp('Example II- with superficial noise');
        TrueImage=SampleImage2;
    end
    
Y = X * TrueImage(:);
noise=randn(ny,1)*Contrast/CNR;
Y=Y+noise;

%The actual model
disp('Computing multiple prior solution');
[lambda,Beta_W,Stats]=nirs_run_DOT_REML(Y,X*W',Beta_prior,Qn,Qp);
%lambda   - hyperparameters
%Beta_W   - the estimated image (in wavelet domain)
%Stats    - model Statistics (in the wavelet domain)

%For comparison, the Tikhonov (equivelent) result
disp('Computing Tikhonov/MNE solution');
[lambda,Beta_W_Tikhonov,Stats]=nirs_run_DOT_REML(Y,X*W',Beta_prior,Qn,QpT);

%Convert to the image domain and display
Beta = W'*Beta_W;
Beta_Tikhonov = W'*Beta_W_Tikhonov;

%Convert the Stats
Stats.tstat.Cov_beta = W'*Stats.tstat.Cov_beta*W;
Stats.tstat.t=Beta./sqrt(diag(Stats.tstat.Cov_beta));
Stats.tstat.pval=2*tcdf(-abs(Stats.tstat.t),Stats.tstat.dfe);


%Now, display the results
Recon_Image=reshape(full(Beta),size(TrueImage));
Recon_Image_Tikhonov=reshape(full(Beta_Tikhonov),size(TrueImage));

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

subplot(3,2,3); hold on;
imagesc(Medium.CompVol.X,Medium.CompVol.Y,squeeze(Recon_Image_Tikhonov(:,:,1,1)));
for idx=1:size(SD.SrcPos,1); text(SD.SrcPos(idx,1),SD.SrcPos(idx,2),['S-' num2str(idx)]); end;
for idx=1:size(SD.DetPos,1); text(SD.DetPos(idx,1),SD.DetPos(idx,2),['D-' num2str(idx)]); end;
caxis([-maxHbO1 maxHbO1]); axis tight; axis off; colorbar;
title('MNE Regularized Layer-I');

subplot(3,2,4); hold on;
imagesc(Medium.CompVol.X,Medium.CompVol.Y,squeeze(Recon_Image_Tikhonov(:,:,2,1)));
for idx=1:size(SD.SrcPos,1); text(SD.SrcPos(idx,1),SD.SrcPos(idx,2),['S-' num2str(idx)]); end;
for idx=1:size(SD.DetPos,1); text(SD.DetPos(idx,1),SD.DetPos(idx,2),['D-' num2str(idx)]); end;
caxis([-maxHbO2 maxHbO2]); axis tight; axis off; colorbar;
title('MNE Regularized Layer-II');

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
end