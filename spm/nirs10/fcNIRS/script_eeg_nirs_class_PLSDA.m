%% script_eeg_nirs_class
%% Choose training and test sets
load('F:\Edgar\Data\NIRS\classifier_data');
nFolds = 10;
total_PPV = zeros([nFolds 1]);
total_sensitivity = zeros([nFolds 1]);
total_specificity = zeros([nFolds 1]);
for iFolds = 1:nFolds
%     dataType = 'EEG+NIRS';
%     dataType = 'NIRS';
%     dataType = 'EEG';
%     dataType = 'HbO';
%     dataType = 'HbR';
%     dataType = 'EEG+HbO';
%     dataType = 'EEG+HbR';
%     dataType = 'Sim';
    
    % Choose 10% of data points as training set
    kFold = 1/nFolds;
    idx_train = unique(randi(nTimePoints, [round(kFold*nTimePoints) 1]));
    % idx_train = 1:round(kFold*nTimePoints);
    % Choose the rest of data points as test set
    idx_test = setdiff(find(1:nTimePoints),idx_train)';
    clc; close all; figure; set (gcf,'color','w')
    set(gcf,'name','Indices')
    subplot(211)
    plot(idx_train,ones(size(idx_train)),'ro');
    hold on; plot(idx_test,ones(size(idx_test)),'bx');
    legend({'Training set' 'Test set'})
    title(sprintf('Training set percentage = %0.1f %%', kFold*100))
    xlabel('Index'); set(gca,'YTick',[])
    % xlim([700 800]);
    ylim([0.98 1.02])
    % Classes (labels) of the training set
    class_train = sz_vector_NIRS_down(idx_train);
    % True classes (labels) of the test set
    class_test = sz_vector_NIRS_down(idx_test);
    
    switch(dataType)
        case 'EEG+NIRS'
            % Training set
            data_train = [dataNIRS_down(idx_train,:) EEGdata_down(idx_train,:)];
            % Testing data set
            data_test = [dataNIRS_down(idx_test,:) EEGdata_down(idx_test,:)];
        case 'NIRS'
            % Training set
            data_train = dataNIRS_down(idx_train,:);
            % Testing data set
            data_test = dataNIRS_down(idx_test,:);
        case 'EEG'
            % Training set
            data_train = EEGdata_down(idx_train,:);
            % Testing data set
            data_test = EEGdata_down(idx_test,:);
        case 'HbO'
            % Training set
            data_train = dataNIRS_down(idx_train,NIRS.Cf.H.C.wl == 1);
            % Testing data set
            data_test = dataNIRS_down(idx_test,NIRS.Cf.H.C.wl == 1);
        case 'HbR'
            % Training set
            data_train = dataNIRS_down(idx_train,NIRS.Cf.H.C.wl == 2);
            % Testing data set
            data_test = dataNIRS_down(idx_test,NIRS.Cf.H.C.wl == 2);
        case 'EEG+HbO'
            % Training set
            data_train = [dataNIRS_down(idx_train,NIRS.Cf.H.C.wl == 1) EEGdata_down(idx_train,:)];
            % Testing data set
            data_test = [dataNIRS_down(idx_test,NIRS.Cf.H.C.wl == 1) EEGdata_down(idx_test,:)];
        case 'EEG+HbR'
            % Training set
            data_train = [dataNIRS_down(idx_train,NIRS.Cf.H.C.wl == 2) EEGdata_down(idx_train,:)];
            % Testing data set
            data_test = [dataNIRS_down(idx_test,NIRS.Cf.H.C.wl == 2) EEGdata_down(idx_test,:)];
        case 'Sim'
            % Set SNR
            SNR = 0.3;
            % Amplitude of the features
            Amp = 1;
            % Number of channles
            nChannelsSim = size(dataNIRS_down,2);
            % nChannelsSim = 10;
            % Simulate features
            simData = Amp*repmat(sz_vector_NIRS_down,[1 nChannelsSim]);
            % Training set
            data_train = simData(idx_train,:);
            % Testing data set
            data_test = simData(idx_test,:);
            
            % Add noise channel by channel (training data)
            for iChannels = 1:size(data_test,2)
                % Computing signal value (RMS)
                signal = sqrt(mean(abs(data_train(:,iChannels)).^2));      % Define signal to be RMS signal
                noiseRatio  = signal / SNR;
                % Noise computation
                noise = noiseRatio * randn(size(data_train(:,iChannels)));
                % Adding noise to training signal
                data_train(:,iChannels) = data_train(:,iChannels) + noise;
            end
            % Add noise channel by channel (data test)
            for iChannels = 1:size(data_train,2)
                % Computing signal value (RMS)
                signal = sqrt(mean(abs(data_test(:,iChannels)).^2));      % Define signal to be RMS signal
                noiseRatio  = signal / SNR;
                % Noise computation
                noise = noiseRatio * randn(size(data_test(:,iChannels)));
                % Adding noise to training signal
                data_test(:,iChannels) = data_test(:,iChannels) + noise;
            end
            
        otherwise
            % Training set
            data_train = [dataNIRS_down(idx_train,:) EEGdata_down(idx_train,:)];
            % Testing data set
            data_test = [dataNIRS_down(idx_test,:) EEGdata_down(idx_test,:)];
    end
    % Plot data sets
    % figure; set (gcf,'color','w')
    % set(gcf,'name','Data sets')
    colormap(flipud(gray(256)))
    subplot(223); imagesc(1:nChannelsSim, idx_train, data_train);
    if strcmp(dataType, 'Sim')
        title(sprintf('Training data (SNR = %0.2f)',SNR));
    else
        title('Training data');
    end
    xlabel('Channels'); ylabel('t(s)')
    ylim([0 nTimePoints])
    subplot(224); imagesc(1:nChannelsSim, idx_test, data_test);
    if strcmp(dataType, 'Sim')
        title(sprintf('Testing data (SNR = %0.2f)',SNR));
    else
        title('Testing data');
    end
    xlabel('Channels'); ylabel('t(s)')
    ylim([0 nTimePoints])
    
    
    %% We calculate the PLSDA model with 2 components by typing:
    tic
    nComps = 2;
    model = plsdafit(data_train,class_train,nComps,'auto','bayes',1);
    % Once the model is calculated, we can see the model performances by typing:
    model.class_param;
    toc
    
    %% Scores, loadings, calculated class, leverages and many other statistics are
    % stored in the model structure. We can proceed by cross validating (with 5
    % venetian blind groups) the PLSDA model with 2 components:
    cv = plsdacv(data_train,class_train,nComps,'auto','vene',5,'bayes');
    disp('Model metrics:')
    disp(cv.class_param)
    
    %% Finally, we can predict the test set samples by using the calibrated model:
    tic
    pred = plsdapred(data_test,model);
    fprintf('Prediction done! \n');
    disp('Prediction metrics:')
    class_param = calc_class_param(pred.class_pred,class_test)
    toc
    
    %% Plot classification results
    
    % Classes
    idx0 = find(pred.class_pred==0);      % Unclassified
    idx1 = find(pred.class_pred==1);      % No seizure
    idx2 = find(pred.class_pred==2);      % Seizure
    % Total negatives
    N_tot = numel(class_test(class_test==1));
    % Total positives
    P_tot = numel(class_test(class_test==2));
    % True negatives
    v1 = zeros(size(class_test));
    v2 = zeros(size(class_test));
    % v1(idx1) = pred.class_pred(idx1);
    % v2(class_test==1) = class_test(class_test==1);
    v1(idx1) = true;
    v2(class_test==1) = true;
    
    TN = sum(v1&v2);
    % True positives
    v3 = zeros(size(class_test));
    v4 = zeros(size(class_test));
    v3(idx2) = pred.class_pred(idx2);
    v4(class_test==2) = class_test(class_test==2);
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
    table2x2 = [TP FP; FN TN];
    fprintf('\tTP: %04d \tFP: %04d\n\tFN: %04d \tTN: %04d\n\n\tPPV: %0.2f\n',table2x2,xCrit.PPV);
    
    % Average values
    total_PPV(iFolds) = xCrit.PPV;
    total_sensitivity(iFolds) = xCrit.sensitivity;
    total_specificity(iFolds) = xCrit.specificity;
    
    % Plot classification
    figure; set(gcf,'color','w')
    % plot(pred.class_pred,'r.')
    hold on
    plot(class_test,'k--')
    ylim([0 2])
    set(gca,'YTick',[0 1 2]);
    set(gca,'YTickLabel',{'Unclassified' 'No Seizure' 'Seizure'});
    xlabel('t (s)','FontSize',12)
    set(gca,'FontSize',12)
    hold on
    % TN
    plot(find(v1&v2), pred.class_pred(v1&v2),'x','Color',[0 127 14]/255)
    % TP
    plot(find(v3&v4), pred.class_pred(v3&v4),'o','Color',[0 127 14]/255)
    % FP
    v5 = zeros(size(class_test));
    v6 = zeros(size(class_test));
    v5(idx1) = true;
    v5(idx0) = true;
    v6(class_test==2) = true;
    plot(find(v5&v6), pred.class_pred(v5&v6), 'ro')
    % FN
    v7 = zeros(size(class_test));
    v8 = zeros(size(class_test));
    v7(idx2) = true;
    v7(idx0) = true;
    v8(class_test==1) = true;
    plot(find(v7&v8), pred.class_pred(v7&v8), 'rx')
    legend({'Real markers' 'True Negatives' 'True Positives' 'False Negatives'  'False Positives'})
end

%% Average values after nFolds
fprintf('\nAverage results for %s: \n\tPPV: %0.2f\n\tsensitivity: %0.2f\n\tspecificity: %0.2f\n',...
    dataType, mean(total_PPV), mean(total_sensitivity), mean(total_specificity));
% EOF

