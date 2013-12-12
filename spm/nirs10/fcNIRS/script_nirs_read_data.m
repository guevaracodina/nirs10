% Script to read NIRS data
%_______________________________________________________________________________
% Copyright (C) 2013 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________
%% Read a random subject
clear; close all; clc
% Le cas 104MAL est extrêmement inhabituel. Il ne pourra pas être utilisé
% pour les analyses.
% NIRSmat = 'F:\Edgar\Data\NIRS\epiNIRS_data\epiNIRS_Processed\epiNIRS\epi104MAL\dataSPMa\coreg\NIRS.mat';
% epi140GG 1st test, nothing special
% NIRSmat = 'F:\Edgar\Data\NIRS\epiNIRS_data\epiNIRS_Processed\epiNIRS\epi140GG\dataSPMa\coreg\NIRS.mat';
% epi127SD Good optodes coverage
NIRSmat = 'F:\Edgar\Data\NIRS\epiNIRS_data\epiNIRS_Processed\epiNIRS\epi127SD\dataSPMa\coreg\NIRS.mat';
% epi143MDM  Good optodes coverage
% NIRSmat = 'F:\Edgar\Data\NIRS\epiNIRS_data\epiNIRS_Processed\epiNIRS\epi143MDM\dataSPMa\coreg\NIRS.mat';
% NIRSmat = 'F:\Edgar\Data\NIRS\epiNIRS_data\ControlNIRS_Processed\ControlNIRS\Sujet001\dataSPMa\coreg\NIRS.mat';
load(NIRSmat)

%% Read all files from the 4th processing level ODtoHbOHbR
% 1) readBOXY
% 2) remove_chn_stdev
% 3) normalize_baseline
% 4) ODtoHbOHbR
preProcStep = 4;
nSessions = numel(NIRS.Dt.fir.pp(preProcStep).p);
dataNIRS = [];
for iFiles = 1:nSessions,
    % currentData [nChannels x nTimePoints]
    currentData = fopen_NIR(NIRS.Dt.fir.pp(preProcStep).p{iFiles}, NIRS.Cf.H.C.N);
    % dataNIRS [nChannels x (nTimePoints x nSessions)]
    dataNIRS = [dataNIRS, currentData];
end
% Time vector
t = linspace(0, (size(dataNIRS,2)-1)/NIRS.Cf.dev.fs , size(dataNIRS,2));

%% Plot data
h = figure; set(h,'color','w')
set(h,'Name',sprintf('Pre-processing step: %s', NIRS.Dt.fir.pp(preProcStep).pre))
subplot(211)
fprintf('Wavelength: %dnm & %dnm\n',NIRS.Cf.dev.wl)
% switch chromophore
%     case 1
%         hb = 'HbO';
%     case 2
%         hb = 'HbR';
%     case 3
%         hb = 'HbT';
% end
hb = 2;
plot(t, dataNIRS(NIRS.Cf.H.C.wl == hb,:)');
ylabel('[\muM]')
xlabel('t [s]')
% title(['\lambda' sprintf(' = %d nm',NIRS.Cf.dev.wl(hb))])
title('HbO')

%% Compute FFT
[fftNIRS, freq] = nirs_positiveFFT(dataNIRS(NIRS.Cf.H.C.wl == hb,:)', NIRS.Cf.dev.fs);
fftNIRSmean = mean(abs(fftNIRS), 2);

%% Plot FFT
figure(h)
subplot(212)
semilogy(freq, fftNIRSmean, 'k-')
xlabel('f [Hz]')
ylabel('Amplitude [a.u.]')

% EOF
