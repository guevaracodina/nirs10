% Script to read NIRS data
%_______________________________________________________________________________
% Copyright (C) 2013 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________
%% Read a random subject
clear; close all; clc
NIRSmat = 'F:\Edgar\Data\NIRS\epiNIRS_data\epiNIRS_Processed\epiNIRS\epi104MAL\dataSPMa\coreg\NIRS.mat';
load(NIRSmat)

%% Read all files from the 4th processing level ODtoHbOHbR
% 1) readBOXY
% 2) remove_chn_stdev
% 3) normalize_baseline
% 4) ODtoHbOHbR
preProcStep = 4;
dataNIRS = [];
for iFiles = 1:numel(NIRS.Dt.fir.pp(preProcStep).p)
    currentData = fopen_NIR(NIRS.Dt.fir.pp(preProcStep).p{iFiles}, NIRS.Cf.H.C.N);
    dataNIRS = [dataNIRS, currentData];
end
% Time vector
t = linspace(0, (size(dataNIRS,2)-1)/NIRS.Cf.dev.fs , size(dataNIRS,2));

%% Plot data
h = figure; set(h,'color','w')
set(h,'Name',sprintf('Pre-processing step: %s', NIRS.Dt.fir.pp(preProcStep).pre))
subplot(211)
fprintf('Wavelength: %dnm & %dnm\n',NIRS.Cf.dev.wl)
% Choose chromophore
%     case 1
%         hb = 'HbO';
%     case 2
%         hb = 'HbR';
%     case 3
%         hb = 'HbT';
% end
hb = 1;
plot(t, dataNIRS(NIRS.Cf.H.C.wl == hb,:)');
ylabel('[\muM]')
xlabel('t [s]')
% title(['\lambda' sprintf(' = %d nm',NIRS.Cf.dev.wl(hb))])
title('HbO ?')

%% Compute FFT
[fftNIRS, freq] = nirs_positiveFFT(dataNIRS(NIRS.Cf.H.C.wl == hb,:)', NIRS.Cf.dev.fs);
fftNIRSmean = mean(abs(fftNIRS), 2);

%% Plot FFT
figure(h)
subplot(212)
loglog(freq, fftNIRSmean)
xlabel('f [Hz]')
ylabel('Amplitude [a.u.]')
% EOF
