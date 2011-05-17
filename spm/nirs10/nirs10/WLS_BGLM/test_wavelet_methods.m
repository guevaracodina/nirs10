function test_wavelet_methods
%test wavelet methods

%get data
d = fopen_NIR('W:\epiNIRSk\epi127SD\dataSPM\hbrepi127SD_004.nir',282);
%get specified SPM model
SPM = [];
load('W:\epiNIRSk\epi127SD\dataSPM\MDLVC282legacy\SPM.mat');
%send channel
Y = d(1,:)';
%choose session
s = 1;
nsess =1;
tSPM = [];
tSPM.Sess = SPM.Sess(s);                
tSPM.xX = SPM.xX;
%find elements of X for session s
nbeta = size(SPM.xX.X,2);
nbetaS = (nbeta-nsess)/nsess;
%last entry is the constant regressor
beta = [(s-1)*nbetaS+1:s*nbetaS nbeta-nsess+s];
svec = SPM.Sess(s).row;                               
K = struct( 'HParam', SPM.xX.K.HParam,...
                'row', svec ,...
                'RT', SPM.xY.RT,...
                'LParam', SPM.xX.K.LParam);
tSPM.xX.K = spm_filter_HPF_LPF_WMDL(K);
tSPM.xX.X = SPM.xX.X(svec,beta);
tSPM.xX.K.row = 1:length(svec);
%for the recent                            
uSPM = tSPM;                            
%for the legacy                           
tSPM = precoloring_batch_legacy(tSPM,Y);
%for the recent
uSPM = precoloring_batch(uSPM,Y);
a=1;
%figure; plot(tSPM.KY); hold on; plot(uSPM.KY,'r'); hold on; plot(tSPM.xX.X(:,1)*50,'g'); hold off
figure; plot(tSPM.xX.pKX(1,:)*50000,'k'); hold on; plot(uSPM.KY,'r'); hold on; 
plot(tSPM.KY,'b'); hold on; plot(tSPM.xX.X(:,1)*50,'g'); hold off
