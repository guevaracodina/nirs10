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
chosenWL = 1;
plot(t, dataNIRS(NIRS.Cf.H.C.wl == chosenWL,:)');
% plot(t, dataNIRS(15,:)');
xlabel('t [s]')
title(['\lambda' sprintf(' = %d nm',NIRS.Cf.dev.wl(chosenWL))])

%% Compute FFT
[fftNIRS, freq] = nirs_positiveFFT(dataNIRS(NIRS.Cf.H.C.wl == chosenWL,:)', NIRS.Cf.dev.fs);
fftNIRSmean = mean(abs(fftNIRS), 2);

%% Plot FFT
figure(h)
subplot(212)
loglog(freq, fftNIRSmean)
xlabel('f [Hz]')
% EOF
