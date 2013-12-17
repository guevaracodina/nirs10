% Script to read NIRS data
%_______________________________________________________________________________
% Copyright (C) 2013 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________
%% Read subject
clear; close all; clc
% epi127SD Good optodes coverage
% NIRSmat = 'F:\Edgar\Data\NIRS\epiNIRS_data\epiNIRS_Processed\epiNIRS\epi127SD\dataSPMa\coreg\NIRS.mat';
% epi143MDM  Good optodes coverage
% NIRSmat = 'F:\Edgar\Data\NIRS\epiNIRS_data\epiNIRS_Processed\epiNIRS\epi143MDM\dataSPMa\coreg\NIRS.mat';
% Control 01
% NIRSmat = 'F:\Edgar\Data\NIRS\epiNIRS_data\ControlNIRS_Processed\ControlNIRS\Sujet001\dataSPMa\coreg\NIRS.mat';
% Control 02
NIRSmat = 'F:\Edgar\Data\NIRS\epiNIRS_data\ControlNIRS_Processed\ControlNIRS\Sujet002\dataSPMa\coreg\NIRS.mat';
load(NIRSmat)


%% Read all files from the 4th processing level ODtoHbOHbR
% 1) readBOXY
% 2) remove_chn_stdev
% 3) normalize_baseline
% 4) ODtoHbOHbR
preProcStep = 4;
nSessions = numel(NIRS.Dt.fir.pp(preProcStep).p);
dataNIRS = [];
for iSessions = 1:nSessions,
    % currentData [nChannels x nTimePoints]
    currentData = fopen_NIR(NIRS.Dt.fir.pp(preProcStep).p{iSessions}, NIRS.Cf.H.C.N);
    % dataNIRS [nChannels x (nTimePoints x nSessions)]
    dataNIRS{iSessions} = currentData;
    % Time vector
    t{iSessions} = linspace(0, (size(currentData,2)-1)/NIRS.Cf.dev.fs , size(currentData,2));
end

%% Filter data
tic
% Band-pass filter configuration
% Filter order
filterOrder = 2;
% Band-pass cut-off frequencies
BPFfreq = [0.009 0.08];
% Filter type
fType = 'butter';
% Passband/Stopband ripple in dB
Rp_Rs = [.5 80];
% Sampling frequency
fs = NIRS.Cf.dev.fs;
[z, p, k] = nirs_temporalBPFconfig(fType, fs, BPFfreq, filterOrder, Rp_Rs);
nChannels = NIRS.Cf.H.C.N;
for iSessions = 1:nSessions,
    for iChannels = 1:nChannels
        dataNIRSfilt{iSessions}(iChannels,:) = temporalBPFrun(squeeze(dataNIRS{iSessions}(iChannels,:)), z, p, k);
    end
end
toc

%% GLM regression
% Global signal only for the time being, despite its issues.
% no short separation measurements available
for iSessions = 1:nSessions,
    brainSignal{iSessions} = mean (dataNIRSfilt{iSessions}, 1);
    % Initialize single voxel 4-D series
    brainSignalRep = zeros([1 1 1 numel(brainSignal)]);
    %% GLM on images here!
    % Constructing inputs required for GLM analysis within the SPM framework
    clear SPM
    V = spm_create_vol(V);
    SPM.xY.VY = spm_vol(PAT.fcPAT.filtNdown.fnameWholeImage{s1, c1});
    y = spm_read_vols(SPM.xY.VY);
    % Preallocating output images
    yRegress = zeros(size(y));
    
    % We take the mean global signal as regressor
    SPM.xX.name = cellstr(['Global Brain Signal']);
    SPM.xX.X = brainSignal{iSessions}';        % Regression is along first dimension. For one regressor it is a column vector.
end


%% Scrubbing

%% Correlation map ROI-to-channels


% EOF
