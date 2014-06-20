%% script_nirs_classification
% Script to read NIRS data and perform PLS classification
%_______________________________________________________________________________
% Copyright (C) 2014 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________
clear; close all; clc;
addpath(genpath('F:\Edgar\Dropbox\PostDoc\NIRS\real_time\classification_toolbox_3.1'))
load('F:\Edgar\Dropbox\PostDoc\NIRS\real_time\NIRS_test_data');
clear k inv_exs2 interpMapHbR interpMapHbO interpMap iSubject iOnsets iFiles Dat...
    EPF2 K K_DC P P_DC Pminus Pminus_DC Q Q_DC Qinterp R R_DC UPDATE_MAP W eTime...
    h1 h2 h3 hLive haxHbO haxHbR hbLim heighthaxHbO heighthaxHbR himHbO himHbR...
    poshaxHbO poshaxHbR screenSize widthhaxHbO widthhaxHbR x_hat x_hat_DC xhat...
    xhat_DC xhatminus xhatminus_DC z eegAmp hb hbLabel;

%% Downsample seizure and NIRS vectors
c = c';
cDec = zeros(round(size(c,1)/round(NIRS.Cf.dev.fs)),size(c,2));
for iChannels = nChannels
    cDec(:,iChannels) = decimate(c(:,iChannels), round(NIRS.Cf.dev.fs));
end
% Categories: No Seizures = 1; Seizures = 2; 
sz_vector = sz_vector + 1;
class_train = round(decimate(sz_vector, round(eeg_fs / (NIRS.Cf.dev.fs/round(NIRS.Cf.dev.fs)))))';
nSamples = numel(class_train);
HbO_train = c(1:nSamples, NIRS.Cf.H.C.wl == 1);   % HbO
% HbR_train = c(1:nSamples, NIRS.Cf.H.C.wl == 2);   % HbR

%% select the number of optimal components by using the plsdacompsel function
HbO_res = plsdacompsel(HbO_train,class_train,'none','vene',5,'bayes');
% HbR_res = plsdacompsel(HbR_train,class_train,'none','vene',5,'bayes');

%% We'll get the error rate in cross validation (and non-error rate in cross
% validation) associated to each component value. Type on the MATLAB command
% window to see the error rates.
HbO_res.er
% HbR_res.er

%% We can then calculate the PLSDA model with 2 components by typing:
tic
model = plsdafit(HbO_train,class_train,2,'none','bayes',1)
toc

%% Predict point-by-point
tic
for i=4:nSamples
    HbO_test = HbO_train(i,:);
    model = plsdafit(HbO_train(1:(i-1),:),class_train(1:(i-1),:),2,'none','bayes',1);
    pred(i) = plsdapred(HbO_test,model);
end
fprintf('Prediction done! \n');
toc
% figure; plot(model.class_calc,'k.')

%% Once the model is calculated, we can see the model performances by typing:
model.class_param

%% Scores, loadings, calculated class, leverages and many other statistics are
% stored in the model structure. We can proceed by cross validating (with 5
% venetian blind groups) the PLSDA model with 2 components:
cv = plsdacv(HbO_train,class_train,2,'none','vene',5,'bayes')

%% Finally, we can predict the test set samples by using the calibrated model:
% Change HbO_train by HbO_test! TO DO...
tic
pred = plsdapred(HbO_train,model)
class_param = calc_class_param(pred.class_pred,class_train)
toc
%% Plot classification results

% Classes
idx0 = find(cv.class_pred==0);      % Unclassified
idx1 = find(cv.class_pred==1);      % No seizure
idx2 = find(cv.class_pred==2);      % Seizure
% Total negatives
N_tot = numel(class_train(class_train==1));
% Total positives
P_tot = numel(class_train(class_train==2));
% True negatives
v1 = zeros(size(class_train));
v2 = zeros(size(class_train));
% v1(idx1) = cv.class_pred(idx1);
% v2(class_train==1) = class_train(class_train==1);
v1(idx1) = true;
v2(class_train==1) = true;

TN = sum(v1&v2);
% True positives
v3 = zeros(size(class_train));
v4 = zeros(size(class_train));
v3(idx2) = cv.class_pred(idx2);
v4(class_train==2) = class_train(class_train==2);
TP = sum(v3&v4);

% False Negatives
FN = P_tot - TP;
% False Positives
FP = N_tot - TN;
% xCrit - Criteria 
xCrit.sensitivity   = TP/(TP + FN); 
xCrit.specificity   = TN/(TN + FP);
% Rate of positive predictions
xCrit.RPP           = (TP+FP)/(TP+FN+FP+TN);
% Rate of negative predictions
xCrit.RNP           = (TN+FN)/(TP+FN+FP+TN);
% Accuracy
xCrit.accu          = (TP+TN)/(TP+FN+FP+TN);
% True positive rate (sensitivity, recall)
xCrit.TPR           = TP/(TP+FN);
% False negative rate (miss)
xCrit.FNR           = FN/(TP+FN);
% False positive rate (fallout)
xCrit.FPR           = FP/(TN+FP);
% True negative rate (specificity)
xCrit.TNR           = TN/(TN+FP);
% Positive predictive value (precision)
xCrit.PPV           = TP/(TP+FP);
% Negative predictive value. 
xCrit.NPV           = TN/(TN+FN);
% Table 2x2
table2x2 = [TP FP; FN TN]

% Plot classification
figure; set(gcf,'color','w')
% plot(cv.class_pred,'r.')
hold on
plot(class_train,'k--')
ylim([0 2])
set(gca,'YTick',[0 1 2]);
set(gca,'YTickLabel',{'Unclassified' 'No Seizure' 'Seizure'});
xlabel('t (s)','FontSize',12)
set(gca,'FontSize',12)
hold on
% TN
plot(find(v1&v2), cv.class_pred(v1&v2),'x','Color',[0 127 14]/255)
% TP
plot(find(v3&v4), cv.class_pred(v3&v4),'o','Color',[0 127 14]/255)
% FP
v5 = zeros(size(class_train));
v6 = zeros(size(class_train));
v5(idx1) = true;
v5(idx0) = true;
v6(class_train==2) = true;
plot(find(v5&v6), cv.class_pred(v5&v6), 'ro')
% FN
v7 = zeros(size(class_train));
v8 = zeros(size(class_train));
v7(idx2) = true;
v7(idx0) = true;
v8(class_train==1) = true;
plot(find(v7&v8), cv.class_pred(v7&v8), 'rx')


legend({'Real' 'True Negatives' 'True Positives' 'False Negatives'  'False Positives'})


%% Kalman figure
% Plot Hb concentrations
h3 = figure;
set(h3, 'color', 'w')
set(h3, 'Name', 'Hemoglobin concentrations')
hold on
% iChannel < 133
iChannel = 98;
plot(t, dataNIRS(iChannel,:), 'r:', 'LineWidth',1);      % HbO
plot(t, dataNIRS(iChannel + nChannels/2,:), 'b:', 'LineWidth',1);      % HbR
plot(t, c(:,iChannel), 'r-', 'LineWidth',2);      % HbO
plot(t, c(:,iChannel + nChannels/2), 'b-', 'LineWidth',2);      % HbR
plot(eeg_t, 700*sz_vector, 'k--', 'LineWidth',2);                % seizure onset
xlim([t(2) t(end)])

title(sprintf('Offline processing'), 'FontSize', 12)
legend({'HbO_2' 'HbR' 'HbO_2 (filt)' 'HbR (filt)' 'Seizure'},...
    'Location', 'SouthEast')
set(gca , 'FontSize', 12)
xlabel('t [s]', 'FontSize', 12); ylabel('Hb [\muM]', 'FontSize', 12); 


% EOF
