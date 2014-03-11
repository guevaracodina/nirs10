%% script_interpolate_channels
clc
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
hb = 2;
dataNIRS = dataNIRS(NIRS.Cf.H.C.wl == hb,:);
% Time vector
t = linspace(0, (size(dataNIRS,2)-1)/NIRS.Cf.dev.fs , size(dataNIRS,2));
% figure; plot(t, dataNIRS');
% Select 1 time point
Dat = dataNIRS(:,100);
figure; plot(t, dataNIRS(50, :));

%% Configuration structure
% Brain view, or a loop can be used here to generate all views
config.brain_view = 2; 
% Option: 0: do not extrapolate
config.AllowExtrapolation = 0; 
% Option: 0: interpolate
config.no_interpolation = 0; 
% Path where interpolated images are saved
config.new_path = 'F:\Edgar\Dropbox\PostDoc\NIRS\real_time'; 
% Name of interpolated image
config.figure_name = 'interp_NIRS_test'; 
% Threshold value of the interpolated image. Attention: Cannot be zero
config.thz = 20; 

%% NIRS data interpolated to cortex
out = nirs_interpolation_render_simplified(Dat, NIRS, config);

% EOF
