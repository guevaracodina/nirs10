%% script_interpolate_channels
clear all; close all; clc
addpath(genpath('D:\spm8\toolbox\nirs10'))
NIRSmat = 'F:\Edgar\Data\NIRS\epiNIRS_data\epiNIRS_Processed\epiNIRS\epi127SD\dataSPMa\coreg\NIRS.mat';
load(NIRSmat)
preProcStep = 4;
% nSessions = numel(NIRS.Dt.fir.pp(preProcStep).p);
dataNIRS = [];
% File number (session)
iFiles = 1;
% currentData [nChannels x nTimePoints]
dataNIRS = fopen_NIR(NIRS.Dt.fir.pp(preProcStep).p{iFiles}, NIRS.Cf.H.C.N);

% Ke's data
% dataNIRS = fopen_NIR('F:\Edgar\Data\NIRS\epiNIRS_data\epiNIRS_Processed\epiNIRS\epi127SD\Ke_Data\epi127SD\Sess1.nir', 3*NIRS.Cf.H.C.N/2);
 
% switch chromophore
%     case 1
%         hb = 'HbO';
%     case 2
%         hb = 'HbR';
%     case 3
%         hb = 'HbT';
% end
hb = 1;
% dataNIRS = dataNIRS(NIRS.Cf.H.C.wl == hb,:);
% Time vector
t = linspace(0, (size(dataNIRS,2)-1)/NIRS.Cf.dev.fs , size(dataNIRS,2));
% figure; plot(t, dataNIRS');
% Select 1 time point
Dat = mean(dataNIRS(:,1000:2000),2);
figure; plot(Dat);

%% Configuration structure
config = nirs_interpolation_render_config (2, 0, 0,...
    'F:\Edgar\Dropbox\PostDoc\NIRS\real_time', 'interp_NIRS_test', 20);

%% NIRS data interpolated to cortex
% Precompute only once
[Q, interpMap, Dat] = nirs_interpolation_render_precompute(Dat, NIRS, config);

% Call nirs_interpolation_render_compute as many times as needed in the loop
tic
interpMap = nirs_interpolation_render_compute(Q, Dat, interpMap);
toc

% old function
% out = nirs_interpolation_render_simplified(Dat, NIRS, config);

%% Display interpolated map

% EOF
