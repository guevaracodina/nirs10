%% script_nirs_classification
% Script to read EEG raw data.
%_______________________________________________________________________________
% Copyright (C) 2014 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________
clear; close all; clc;
% 1) List of all NIRS matrices
NIRSmat{1} = 'epi101LH';
NIRSmat{2} = 'epi102FA';
NIRSmat{3} = 'epi103GJ';
% Le cas 104MAL est extrêmement inhabituel. Il ne pourra pas être utilisé
% pour les analyses. (M)
NIRSmat{4} = 'F:\Edgar\Data\NIRS\epiNIRS_data\epiNIRS_Processed\epiNIRS\epi104MAL\dataSPMa\coreg\NIRS.mat';
NIRSmat{5} = 'epi105MER';
NIRSmat{6} = 'epi106VLL';
NIRSmat{7} = 'epi107GC';
NIRSmat{8} = 'epi108DD';
NIRSmat{9} = 'epi109OC';
NIRSmat{10} = 'epiLFCDB'; %PP put this subject here for convenience
NIRSmat{11} = 'epi111ML';
NIRSmat{12} = 'epi112VV';
NIRSmat{13} = 'epi113AG';
NIRSmat{14} = 'epi114MCG';
NIRSmat{15} = 'epi115MD';
NIRSmat{16} = 'epi116GA';
NIRSmat{17} = 'epi117FDH';
NIRSmat{18} = 'epi118MC';
NIRSmat{19} = 'epi119MG';
NIRSmat{21} = 'epi121FB';
% epi122 - Seizures in sessions 1,3,4,6
NIRSmat{22} = 'F:\Edgar\Data\NIRS\epiNIRS_data\epiNIRS_Processed\epiNIRS\epi122PEV\dataSPMa\coreg\NIRS.mat';
NIRSmat{23} = 'epi123JR';
NIRSmat{24} = 'epi124YP';
NIRSmat{25} = 'epi125YD';
NIRSmat{26} = 'epi126AB';
% epi127SD Good optodes coverage
NIRSmat{27} = 'F:\Edgar\Data\NIRS\epiNIRS_data\epiNIRS_Processed\epiNIRS\epi127SD\dataSPMa\coreg\NIRS.mat';
NIRSmat{28} = 'epi128RW';
NIRSmat{29} = 'epi129DP';
NIRSmat{30} = 'epi130CO';
NIRSmat{31} = 'epi131AMA';
NIRSmat{32} = 'epi132FG';
NIRSmat{33} = 'epi133EA';
NIRSmat{34} = 'epi134ASC';
NIRSmat{35} = 'epi135MHT';
NIRSmat{36} = 'epi136JSL';
NIRSmat{37} = 'epi137SR';
NIRSmat{38} = 'epi138BG';
NIRSmat{39} = 'epi139JB';
% epi140GG 1st test, nothing special
NIRSmat{40} = 'F:\Edgar\Data\NIRS\epiNIRS_data\epiNIRS_Processed\epiNIRS\epi140GG\dataSPMa\coreg\NIRS.mat';
NIRSmat{41} = 'epi141JP';
NIRSmat{42} = 'epi142EG';
% epi143MDM  Good optodes coverage
NIRSmat{43} = 'F:\Edgar\Data\NIRS\epiNIRS_data\epiNIRS_Processed\epiNIRS\epi143MDM\dataSPMa\coreg\NIRS.mat';
NIRSmat{44} = 'epi144TL';
NIRSmat{45} = 'epi145FS';
NIRSmat{46} = 'epi146MM'; %No T1 available, used epi127SD to run script, but only project on template
NIRSmat{47} = 'epi147LL';
NIRSmat{48} = 'epi148MHC';
NIRSmat{49} = 'epi149AHL';
nSubjects = length(NIRSmat);
% Choose subject to analyze
iSubject = 22;
% Load info
load(NIRSmat{iSubject})
addpath(genpath('F:\Edgar\Dropbox\Matlab'))
addpath(genpath('F:\Edgar\Dropbox\PostDoc\NIRS\real_time\classification_toolbox_3.1'))
% load('F:\Edgar\Dropbox\PostDoc\NIRS\real_time\NIRS_test_data');
% clear k inv_exs2 interpMapHbR interpMapHbO interpMap iSubject iOnsets iFiles Dat...
%     EPF2 K K_DC P P_DC Pminus Pminus_DC Q Q_DC Qinterp R R_DC UPDATE_MAP W eTime...
%     h1 h2 h3 hLive haxHbO haxHbR hbLim heighthaxHbO heighthaxHbR himHbO himHbR...
%     poshaxHbO poshaxHbR screenSize widthhaxHbO widthhaxHbR x_hat x_hat_DC xhat...
%     xhat_DC xhatminus xhatminus_DC z eegAmp hb hbLabel;

%% Read all files from the 4th processing level ODtoHbOHbR
% 1) readBOXY
% 2) remove_chn_stdev
% 3) normalize_baseline
% 4) ODtoHbOHbR
preProcStep = 4;
nSessions = numel(NIRS.Dt.fir.pp(preProcStep).p);
% Choose file(session)
% iFiles = 3;
% Preallocation
dataNIRS = cell([1 nSubjects]);
t = cell([1 nSubjects]);
for iFiles = 1:numel(NIRS.Dt.fir.Sess)
    fileName = NIRS.Dt.fir.pp(preProcStep).p{iFiles};
    fprintf('Reading %s...\n', fileName);
    if preProcStep == 4
        dataNIRS{iSubject}{iFiles} = fopen_NIR(fileName, NIRS.Cf.H.C.N)';
    else
        dataNIRS{iSubject}{iFiles} = fopen_NIR(fileName, NIRS.Cf.H.C.id(1,end))';
    end
    % Time vector
    t{iSubject}{iFiles} = linspace(0, (size(dataNIRS{iSubject}{iFiles},1)-1)/NIRS.Cf.dev.fs , size(dataNIRS{iSubject}{iFiles},1))';
    fprintf('Reading %s done!\n', fileName);
end

%% Open ASCII .dat file
pathName = 'F:\Edgar\Data\NIRS\epiNIRS_data\epiNIRS\epi122PEV\eeg\raw';
fName{iSubject}{1} = 'epi122_C1_Raw Data';
fName{iSubject}{2} = 'epi122_C2_Raw Data';
fName{iSubject}{3} = 'epi122_C3_Raw Data';
fName{iSubject}{4} = 'epi122_C4_Raw Data';
fName{iSubject}{5} = 'epi122_C5_Raw Data';
fName{iSubject}{6} = 'epi122_C6_Raw Data';
% Number of data points
dataPoints{iSubject}{1} = 465800;
dataPoints{iSubject}{2} = 458480;
dataPoints{iSubject}{3} = 456100;
dataPoints{iSubject}{4} = 456540;
dataPoints{iSubject}{5} = 464120;
dataPoints{iSubject}{6} = 215440;
% Number of channels
nChannelsEEG{iSubject}{1} = 22;
nChannelsEEG{iSubject}{2} = 22;
nChannelsEEG{iSubject}{3} = 22;
nChannelsEEG{iSubject}{4} = 22;
nChannelsEEG{iSubject}{5} = 22;
nChannelsEEG{iSubject}{6} = 22;
% Sampling interval (in us)
SamplingInterval = 2000;
% Sampling interval (in us)
TR = SamplingInterval/1e6;
% EEG sampling frequency (Hz)
eeg_fs = 1 / TR;
% Preallocation
EEGdata = cell([1 nSubjects]);
EEGdata{iSubject} = cell([1 length(fName{iSubject})]);
EEGtmp = cell([size(fName,1) nChannelsEEG{iSubject}{1}]);
hdr = cell([size(fName,1) nChannelsEEG{iSubject}{1}]);
for iFiles = 1:length(fName{iSubject})
    fileName = fullfile(pathName,[fName{iSubject}{iFiles} '.dat']);
    fid = fopen(fileName);
    fprintf('Reading %s...\n', fileName);
    for iChannels = 1:nChannelsEEG{iSubject}{iFiles}
        hdr(iFiles, iChannels) = textscan(fid,'%s', 1);
        EEGtmp(iFiles, iChannels) = textscan(fid, '%f', dataPoints{iSubject}{iFiles}, 'CollectOutput', 1);
        EEGdata{iSubject}{iFiles} = [EEGdata{iSubject}{iFiles} EEGtmp{iFiles, iChannels}];
    end
    % get rid of VeoG, EKG and Remg (i.e. take only first 19 channels)
    EEGdata{iSubject}{iFiles} = EEGdata{iSubject}{iFiles}(:,1:19);
    fprintf('Reading %s done!\n', fileName);
end
fclose(fid);
EEGhdr{iSubject} = hdr(1,1:19);
% cleanup
clear EEGtmp hdr

%% Seizure markers (labels) vector
% Epi122PEV (Session 1,3,4,6) F:\Edgar\Data\NIRS\epiNIRS_data\epiNIRS\epi122PEV\eeg\szMarkers
% \epi122_C1_AllMarkers (25 hits)
% Line 63: Stimulus, sz start, 19486, 1, F7
% Line 64: Stimulus, szcontinu, 19503, 2491, T3
% Line 83: Stimulus, szcontinu, 21997, 4497, F7
% Line 118: Stimulus, szcontinu, 26420, 4574, T5
% Line 156: Stimulus, szcontinu, 31002, 4492, F3
% Line 188: Stimulus, szcontinu, 35517, 4477, P3
% Line 201: Stimulus, szcontinu, 39997, 4497, F3
% Line 206: Stimulus, szcontinu, 44492, 4502, P3
% Line 208: Stimulus, szcontinu, 49007, 4487, F3
% Line 209: Stimulus, szcontinu, 53456, 4538, P3
% Line 214: Stimulus, szcontinu, 57987, 4221, C3
% Line 226: Stimulus, szcontinu, 62215, 814, P3
% Line 229: Stimulus, sz end, 63000, 1, F7
% Line 1218: Stimulus, szcontinu, 419416, 2578, T3
% Line 1219: Stimulus, sz start, 419437, 1, F7
% Line 1226: Stimulus, szcontinu, 421976, 4518, T5
% Line 1239: Stimulus, szcontinu, 426517, 4477, T3
% Line 1254: Stimulus, szcontinu, 430981, 4513, T5
% Line 1270: Stimulus, szcontinu, 435461, 4533, F3
% Line 1289: Stimulus, szcontinu, 439976, 4518, P3
% Line 1310: Stimulus, szcontinu, 444471, 4523, F3
% Line 1328: Stimulus, szcontinu, 449033, 4461, Fp1
% Line 1343: Stimulus, szcontinu, 453461, 4533, T5
% Line 1358: Stimulus, szcontinu, 457987, 4507, F3
% Line 1373: Stimulus, szcontinu, 462466, 3331, C3
sz_startIdx{iSubject}{1}    = [19486 419437];
sz_endIdx{iSubject}{1}      = [63000 462466+3331];
% \epi122_C3_AllMarkers (15 hits)
% Line 389: Stimulus, szcontinu, 107669, 2326, T5
% Line 390: Stimulus, sz start, 107704, 1, F7
% Line 398: Stimulus, szcontinu, 109992, 4503, Fp1
% Line 412: Stimulus, szcontinu, 114497, 3792, Fp1
% Line 427: Stimulus, szcontinu, 118328, 4667, T5
% Line 446: Stimulus, szcontinu, 122973, 4522, P3
% Line 466: Stimulus, szcontinu, 127502, 4493, F3
% Line 487: Stimulus, szcontinu, 131987, 4508, Fp1
% Line 510: Stimulus, szcontinu, 136487, 4994, P3
% Line 534: Stimulus, szcontinu, 141455, 3980, F3
% Line 552: Stimulus, szcontinu, 145403, 4592, C3
% Line 570: Stimulus, szcontinu, 149953, 4542, Fp2
% Line 587: Stimulus, szcontinu, 154492, 4503, P3
% Line 604: Stimulus, szcontinu, 158987, 1312, C3
% Line 609: Stimulus, sz end, 160302, 1, F7
sz_startIdx{iSubject}{3}    = 107704;
sz_endIdx{iSubject}{3}      = 160302;
% \epi122_C4_AllMarkers (11 hits)
% Line 743: Stimulus, sz start, 229124, 1, F7
% Line 744: Stimulus, szcontinu, 229177, 4317, T5
% Line 759: Stimulus, szcontinu, 233466, 4528, T3
% Line 776: Stimulus, szcontinu, 238040, 4454, T5
% Line 777: Stimulus, szcontinu, 238076, 1, F7
% Line 795: Stimulus, szcontinu, 242370, 4624, F7
% Line 815: Stimulus, szcontinu, 247013, 4481, T5
% Line 833: Stimulus, szcontinu, 251286, 4708, F3
% Line 850: Stimulus, szcontinu, 255981, 4513, T5
% Line 865: Stimulus, szcontinu, 260503, 2721, Fp1
% Line 874: Stimulus, sz end, 263207, 1, F7
sz_startIdx{iSubject}{4}    = 229124;
sz_endIdx{iSubject}{4}      = 263207;
% \epi122_C6_AllMarkers (13 hits)
% Line 337: Stimulus, sz start, 95117, 1, F7
% Line 338: Stimulus, szcontinu, 95134, 4360, T3
% Line 351: Stimulus, szcontinu, 99496, 4498, T5
% Line 366: Stimulus, szcontinu, 103991, 4503, F7
% Line 381: Stimulus, szcontinu, 108266, 4728, O2
% Line 399: Stimulus, szcontinu, 112909, 4585, T3
% Line 418: Stimulus, szcontinu, 117240, 4754, F4
% Line 441: Stimulus, szcontinu, 121807, 4687, T3
% Line 464: Stimulus, szcontinu, 126246, 4748, F4
% Line 486: Stimulus, szcontinu, 130710, 4784, T3
% Line 508: Stimulus, szcontinu, 135486, 4508, T6
% Line 527: Stimulus, szcontinu, 139899, 577, F4
% Line 530: Stimulus, sz end, 140462, 1, F7
sz_startIdx{iSubject}{6}    = 95117;
sz_endIdx{iSubject}{6}      = 140462;

for iFiles = 1:length(fName{iSubject})
% iFiles=1
    % time vector for eeg
%     eeg_t{iSubject}{iFiles} = linspace(t{iSubject}{iFiles}(1), t{iSubject}{iFiles}(end), round(eeg_fs*(t{iSubject}{iFiles}(end) - t{iSubject}{iFiles}(1))));
    eeg_t{iSubject}{iFiles} = linspace(t{iSubject}{iFiles}(1), t{iSubject}{iFiles}(end), size(EEGdata{iSubject}{iFiles},1))';
    % Seizure vector
    sz_vector{iSubject}{iFiles} = zeros(size(eeg_t{iSubject}{iFiles}));
    sz_vector_NIRS{iSubject}{iFiles} = zeros(size(t{iSubject}{iFiles}));
    % Seizure start index
    sz_start{iSubject}{iFiles} = sz_startIdx{iSubject}{iFiles};
    % Seizure end index
    sz_end{iSubject}{iFiles} = sz_endIdx{iSubject}{iFiles};
    for iOnsets = 1:numel(sz_start{iSubject}{iFiles})
        sz_vector{iSubject}{iFiles}(sz_start{iSubject}{iFiles}(iOnsets):sz_end{iSubject}{iFiles}(iOnsets)) = 1;
        sz_vector_NIRS{iSubject}{iFiles} = downsample(sz_vector{iSubject}{iFiles}, round(eeg_fs/NIRS.Cf.dev.fs));
        if size(sz_vector_NIRS{iSubject}{iFiles},1) >= size(dataNIRS{iSubject}{iFiles},1)
            sz_vector_NIRS{iSubject}{iFiles} = sz_vector_NIRS{iSubject}{iFiles}(1:size(dataNIRS{iSubject}{iFiles},1),:);
        else
            % append zeros at the end
            zeros2append = zeros([(size(dataNIRS{iSubject}{iFiles},1)-size(sz_vector_NIRS{iSubject}{iFiles},1)), 1]);
            sz_vector_NIRS{iSubject}{iFiles} = [sz_vector_NIRS{iSubject}{iFiles}; zeros2append];
        end
    end
end
% Amplitud of seizure onset vector, to be displayes correctly
eegAmp = 1000;
% figure; stem(eeg_t, sz_vector)

%% Plot EEG and seizures for EPI122
yLimits = [-1000 1000];
figure; set(gcf,'color','w')
% iFiles = 6;
% The size of the matrix plot
M = 3;
N = 4;
for iFiles = 1:length(fName{iSubject})
    % The index according to your preferred ordering (column-wise)
    % i_colwise = 4;
    % Conversion function
    [jj,ii] = ind2sub([N,M],2*iFiles-1);
    i_rowwise = sub2ind([M,N],ii,jj); % This is the ordering MATLAB expects (row-wise)
    subplot (N,M,i_rowwise)
    plot(eeg_t{iSubject}{iFiles},EEGdata{iSubject}{iFiles},'k-')
    ylim(yLimits)
    ylabel('EEG (\muV)','FontSize',12)
    xlabel('t (s)','FontSize',12)
    [jj,ii] = ind2sub([N,M],2*iFiles);
    i_rowwise = sub2ind([M,N],ii,jj); % This is the ordering MATLAB expects (row-wise)
    subplot(N,M,i_rowwise)
    plot(eeg_t{iSubject}{iFiles},sz_vector{iSubject}{iFiles},'k:','LineWidth',2)
    ylim([-0.1 1.1])
    xlabel('t (s)','FontSize',12)
    ylabel('Seizure','FontSize',12)
    set(gca,'YTick',[0 1]); set(gca,'YTickLabel',{'No' 'Yes'})
end

%% Plot NIRS and seizures for EPI122
yLimits = [-1000 1000];
figure; set(gcf,'color','w')
% iFiles = 6;
% The size of the matrix plot
M = 3;
N = 4;
for iFiles = 1:length(fName{iSubject})
    % The index according to your preferred ordering (column-wise)
    % i_colwise = 4;
    % Conversion function
    [jj,ii] = ind2sub([N,M],2*iFiles-1);
    i_rowwise = sub2ind([M,N],ii,jj); % This is the ordering MATLAB expects (row-wise)
    subplot (N,M,i_rowwise)
    hold on
    % HbO=1
    HbOdata = dataNIRS{iSubject}{iFiles}(:, NIRS.Cf.H.C.wl == 1);
    % HbR=2
    HbRdata = dataNIRS{iSubject}{iFiles}(:, NIRS.Cf.H.C.wl == 2);
    plot(t{iSubject}{iFiles},HbOdata,'r-')
    plot(t{iSubject}{iFiles},HbRdata,'b-')
    ylim(yLimits)
    ylabel('\DeltaHb (\muM)','FontSize',12)
    xlabel('t (s)','FontSize',12)
    [jj,ii] = ind2sub([N,M],2*iFiles);
    i_rowwise = sub2ind([M,N],ii,jj); % This is the ordering MATLAB expects (row-wise)
    subplot(N,M,i_rowwise)
    plot(t{iSubject}{iFiles},sz_vector_NIRS{iSubject}{iFiles},'k:','LineWidth',2)
    ylim([-0.1 1.1])
    xlabel('t (s)','FontSize',12)
    ylabel('Seizure','FontSize',12)
    set(gca,'YTick',[0 1]); set(gca,'YTickLabel',{'No' 'Yes'})
end

%% Joining all data
dataNIRS_total = [];
sz_vector_NIRS_total = [];
EEGdata_total = [];
sz_vector_total = [];
for iFiles = 1:length(fName{iSubject})
    dataNIRS_total = [dataNIRS_total; dataNIRS{iSubject}{iFiles}];
    sz_vector_NIRS_total = [sz_vector_NIRS_total; sz_vector_NIRS{iSubject}{iFiles}];
    EEGdata_total = [EEGdata_total; EEGdata{iSubject}{iFiles}];
    sz_vector_total = [sz_vector_total; sz_vector{iSubject}{iFiles}];
end

%% Plot amplitude distributions
% HbO
figure; set(gcf,'color','w')
subplot(221)
nhist({dataNIRS_total(sz_vector_NIRS_total==0,NIRS.Cf.H.C.wl == 1), dataNIRS_total(sz_vector_NIRS_total==1,NIRS.Cf.H.C.wl == 1)},...
    'pdf','boxplot','legend',{'No Sz', 'Sz'},'color','sequential')
title({'HbO_2'})
ylim([0 0.04 ]);
xlim([-500 1000]);
subplot(222)
nhist({dataNIRS_total(sz_vector_NIRS_total==0,NIRS.Cf.H.C.wl == 2), dataNIRS_total(sz_vector_NIRS_total==1,NIRS.Cf.H.C.wl == 2)},...
    'pdf','boxplot','legend',{'No Sz', 'Sz'},'color','sequential')
title({'HbR'})
ylim([0 0.04 ]);
xlim([-500 1000]);
subplot(212)
nhist({EEGdata_total(sz_vector_total==0,:), EEGdata_total(sz_vector_total==1,:)},...
    'pdf','boxplot','legend',{'No Sz', 'Sz'},'color','sequential')
title({'EEG'})


%% Downsample seizure and NIRS vectors
% dataNIRS_down = zeros(round(size(dataNIRS_total,1)/NIRS.Cf.dev.fs),size(dataNIRS_total,2));
% sz_vector_NIRS_down = zeros(round(size(sz_vector_NIRS_total,1)/NIRS.Cf.dev.fs),size(sz_vector_NIRS_total,2));
% EEGdata_down = zeros(round(size(EEGdata_total,1)/eeg_fs),size(EEGdata_total,2));
% sz_vector_down = zeros(round(size(sz_vector_total,1)/eeg_fs),size(sz_vector_total,2));
clear dataNIRS_down sz_vector_NIRS_down EEGdata_down EEGdata_down
for iChannels=1:size(dataNIRS_total,2)
    dataNIRS_down(:,iChannels) = decimate(dataNIRS_total(:,iChannels), round(NIRS.Cf.dev.fs));
end
sz_vector_NIRS_down = round(decimate(sz_vector_NIRS_total(:,1), round(NIRS.Cf.dev.fs)));
for iChannels=1:size(EEGdata_total,2)
    EEGdata_down(:,iChannels) = decimate(EEGdata_total(:,iChannels), round(eeg_fs));
end
sz_vector_down = round(decimate(sz_vector_total(:,1), round(eeg_fs)));

% Find minimum number of time points to keep
nTimePoints = min([size(sz_vector_down,1) size(EEGdata_down,1) size(sz_vector_NIRS_down,1) size(dataNIRS_down,1)]);
dataNIRS_down = dataNIRS_down(1:nTimePoints,:);
sz_vector_NIRS_down = sz_vector_NIRS_down(1:nTimePoints,:);
EEGdata_down = EEGdata_down(1:nTimePoints,:);
sz_vector_down = sz_vector_down(1:nTimePoints,:);
% Categories: No Seizures = 1; Seizures = 2; 
sz_vector_NIRS_down = sz_vector_NIRS_down + 1;
figure; stem(sz_vector_NIRS_down)

%%
% c = c';
% cDec = zeros(round(size(c,1)/round(NIRS.Cf.dev.fs)),size(c,2));
% for iChannels = nChannels{iSubject}
%     cDec(:,iChannels) = decimate(c(:,iChannels), round(NIRS.Cf.dev.fs));
% end

% sz_vector = sz_vector + 1;
% class_train = round(decimate(sz_vector, round(eeg_fs / (NIRS.Cf.dev.fs/round(NIRS.Cf.dev.fs)))))';
% nSamples = numel(class_train);
% HbO_train = c(1:nSamples, NIRS.Cf.H.C.wl == 1);   % HbO
% HbR_train = c(1:nSamples, NIRS.Cf.H.C.wl == 2);   % HbR

%% Choose training and test sets
% Choose 10% of data points as training set
idx_train = unique(randi(nTimePoints, [round(nTimePoints/10) 1]));
% Choose the rest of data points as test set
idx_test = setdiff(find(1:nTimePoints),idx_train);
class_train = sz_vector_NIRS_down(idx_train);
data_train = [dataNIRS_down(idx_train,:) EEGdata_down(idx_train,:)];
class_test = sz_vector_NIRS_down(idx_test);
data_test = [dataNIRS_down(idx_test,:) EEGdata_down(idx_test,:)];

%% select the number of optimal components by using the plsdacompsel function
data_res = plsdacompsel(data_train,class_train,'none','vene',5,'bayes');
% HbR_res = plsdacompsel(HbR_train,class_train,'none','vene',5,'bayes');

%% We'll get the error rate in cross validation (and non-error rate in cross
% validation) associated to each component value. Type on the MATLAB command
% window to see the error rates.
data_res.er
% HbR_res.er

%% We can then calculate the PLSDA model with 2 components by typing:
tic
model = plsdafit(data_train,class_train,2,'none','bayes',1)
toc

%% Predict point-by-point (Do ten 10-fold cross validation run)
% tic
% for i=4:nSamples
%     HbO_test = HbO_train(i,:);
%     model = plsdafit(HbO_train(1:(i-1),:),class_train(1:(i-1),:),2,'none','bayes',1);
%     pred = plsdapred(data_test,model);
% end

% toc
% figure; plot(model.class_calc,'k.')

%% Once the model is calculated, we can see the model performances by typing:
model.class_param

%% Scores, loadings, calculated class, leverages and many other statistics are
% stored in the model structure. We can proceed by cross validating (with 5
% venetian blind groups) the PLSDA model with 2 components:
cv = plsdacv(data_train,class_train,2,'none','vene',5,'bayes')

%% Finally, we can predict the test set samples by using the calibrated model:
% Change HbO_train by HbO_test! TO DO...
tic
pred = plsdapred(data_test,model)
class_param = calc_class_param(pred.class_pred,class_test)
fprintf('Prediction done! \n');
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
% % Plot Hb concentrations
% h3 = figure;
% set(h3, 'color', 'w')
% set(h3, 'Name', 'Hemoglobin concentrations')
% hold on
% % iChannel < 133
% iChannel = 98;
% plot(t, dataNIRS(iChannel,:), 'r:', 'LineWidth',1);      % HbO
% plot(t, dataNIRS(iChannel + nChannels/2,:), 'b:', 'LineWidth',1);      % HbR
% plot(t, c(:,iChannel), 'r-', 'LineWidth',2);      % HbO
% plot(t, c(:,iChannel + nChannels/2), 'b-', 'LineWidth',2);      % HbR
% plot(eeg_t, 700*sz_vector, 'k--', 'LineWidth',2);                % seizure onset
% xlim([t(2) t(end)])
% 
% title(sprintf('Offline processing'), 'FontSize', 12)
% legend({'HbO_2' 'HbR' 'HbO_2 (filt)' 'HbR (filt)' 'Seizure'},...
%     'Location', 'SouthEast')
% set(gca , 'FontSize', 12)
% xlabel('t [s]', 'FontSize', 12); ylabel('Hb [\muM]', 'FontSize', 12); 


% EOF

