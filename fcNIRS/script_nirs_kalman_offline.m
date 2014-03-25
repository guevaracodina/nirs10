% Script to read NIRS data and perform Kalman filtering offline
%_______________________________________________________________________________
% Copyright (C) 2014 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________
%% Read a subject
% clear; close all; clc
% NIRSmat = 'F:\Edgar\Data\NIRS\epiNIRS_data\ControlNIRS_Processed\ControlNIRS\Sujet001\dataSPMa\coreg\NIRS.mat';

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
% Choose subject to analyze
iSubject = 22;
% Load info
load(NIRSmat{iSubject})

%% Read all files from the 4th processing level ODtoHbOHbR
% 1) readBOXY
% 2) remove_chn_stdev
% 3) normalize_baseline
% 4) ODtoHbOHbR
preProcStep = 1;
nSessions = numel(NIRS.Dt.fir.pp(preProcStep).p);
% Choose file(session)
iFiles = 3;
if preProcStep == 4
    dataNIRS = fopen_NIR(NIRS.Dt.fir.pp(preProcStep).p{iFiles}, NIRS.Cf.H.C.N);
else
    dataNIRS = fopen_NIR(NIRS.Dt.fir.pp(preProcStep).p{iFiles},NIRS.Cf.H.C.id(1,end));
end
% Time vector
t = linspace(0, (size(dataNIRS,2)-1)/NIRS.Cf.dev.fs , size(dataNIRS,2));

%% EEG seizures list
% Epi122PEV (Session 1,3,4,6)
sz_startIdx{22}{3} = 107704;
sz_endIdx{22}{3} = 160302;

% eeg sampling frequency
eeg_fs = 1/2e-3;
% time vector for eeg
eeg_t = linspace(t(1), t(end), round(eeg_fs*(t(end) - t(1))));
% Seizure vector
sz_vector = zeros(size(eeg_t));
% Seizure start index
sz_start = sz_startIdx{iSubject}{iFiles};
% Seizure end index
sz_end = sz_endIdx{iSubject}{iFiles};
for iOnsets = 1:numel(sz_start)
    sz_vector(sz_start(iOnsets):sz_end(iOnsets)) = 1;
end
% Amplitud of seizure onset vector, to be displayes correctly
eegAmp = 1000;
% figure; stem(eeg_t, sz_vector)

%% Plot raw data
h1 = figure; set(h1,'color','w')
set(h1,'Name',sprintf('Pre-processing step: %s', NIRS.Dt.fir.pp(preProcStep).pre))
subplot(211)
fprintf('Wavelength: %dnm & %dnm\n',NIRS.Cf.dev.wl)
hb = 2;
switch hb
    case 1
        hbLabel = 'HbO';
    case 2
        hbLabel = 'HbR';
    case 3
        hbLabel = 'HbT';
end
plot(t, dataNIRS(NIRS.Cf.H.C.wl == hb,:)');
ylabel('Raw amplitude [a.u.]')
xlabel('t [s]')

%% Kalman filtering initial settings
% Numer of iterations ()
nIter = numel(t);
% Measurement (Observations)
% z = mean(dataNIRS(NIRS.Cf.H.C.wl == hb,:));     % Mean of channels
% z = dataNIRS(NIRS.Cf.H.C.wl == hb,:);           % Channels at 1 wavelength
z = dataNIRS;                                   % All Channels
z = z(1:NIRS.Cf.H.C.N,:);                       % Except the last 2 channels
sizeMat = size(z);
nChannels = sizeMat(1);
%  Very small process variance:
Q = 1e-5*ones([nChannels 1]);
% allocate space for arrays
xhat = zeros(sizeMat);          % a posteriori estimate of x
P = zeros(sizeMat);             % a posteriori error estimate
xhatminus = zeros(sizeMat);     % a priori estimate of x
Pminus = zeros(sizeMat);        % a priori error estimate
K = zeros(sizeMat);             % gain or blending factor

% Initial guess the measurement is zero
x_hat(:,1) = 0;
%  Initial guess Choose P_0 to be 1
P(:,1) = 1;

% R is the estimate of measurement variance, change to see effect
R = 0.01*ones([nChannels 1]);
% R = 2.2*ones([nChannels 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DC value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Very small process variance:
Q_DC = 1e-5*ones([nChannels 1]);
% allocate space for arrays
xhat_DC = zeros(sizeMat);           % a posteriori estimate of x
P_DC = zeros(sizeMat);              % a posteriori error estimate
xhatminus_DC = zeros(sizeMat);      % a priori estimate of x
Pminus_DC = zeros(sizeMat);         % a priori error estimate
K_DC = zeros(sizeMat);              % gain or blending factor
OD = zeros(sizeMat);                % Optical Density
c = zeros(sizeMat);                 % Concentrations
% Initial guess the measurement is zero
x_hat_DC(:,1) = 0;
%  Initial guess Choose P_0 to be 1
P_DC(:,1) = 1;

% R is the estimate of measurement variance, change to see effect
R_DC = 10*ones([nChannels 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hemoglobin concentrations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partial volume correction factor
PVF = [50 50];
% Differential pathlength factor
% From Duncan et al., Phys. Med. Biol. 1995 40(2):295-304.
DPF = [6.51 5.86]; % For 690nm and 832nm
% NOTE: Portable NIRS-EEG wavelengths are 735nm and 850nm
DPF = repmat(DPF',[1 nChannels]); % nLambda x nChannels
PVF = repmat(PVF',[1 nChannels]); % nLambda x nChannels
EPF = zeros(1, nChannels); %EPF = L * PPF = L * DPF/PVF at each wavelength
% Source-detector distances
Cgp = NIRS.Cf.H.C.gp;
% Channels wavelength
Cwl = NIRS.Cf.H.C.wl;
% Device wavelength
wl = NIRS.Cf.dev.wl;
% for iChannels = 1:nChannels
for iChannels = 1:NIRS.Cf.H.C.N
    EPF(1,iChannels) = Cgp(1,iChannels)*DPF(Cwl(1,iChannels),iChannels)./PVF(Cwl(1,iChannels),iChannels);
end
% Effective path length
EPF2 = ones(nIter,1)*EPF; %PP used to be EPF2 = EPF *ones(1,size(d,2));
% exs(:,1): HbO for each wavelength; exs(:,2): HbR for each wavelength Alexis's
% choice of extinction coefficients corresponds to case 1 in Homer, which
% appears to be their preferred choice too
[exs, extcoeff_ref] = GetExtinctions(wl,1);
inv_exs = pinv(exs(:,1:2)); % inv_exs has size 2 x #wl
inv_exs2 = kron(inv_exs,eye(nChannels/size(wl,2))); % size #pairs*2 x #pairs*#wl
% (number of channels nChannels = number of pairs x number of wavelengths)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolation onto the cortex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration structure
config = nirs_interpolation_render_config (2, 0, 0,...
    'F:\Edgar\Dropbox\PostDoc\NIRS\real_time', 'interp_NIRS_test', 20);
% NIRS data interpolated to cortex
% Precompute only once for each contrast
[Qinterp, interpMap, Dat, W] = nirs_interpolation_render_precompute(z(NIRS.Cf.H.C.wl == 1, 1), NIRS, config);
interpMapHbR = interpMap;
interpMapHbO = interpMap;
% Preallocate interpMapMovie (very big matrix)
% interpMapMovieHbO = zeros([size(interpMap) nIter]);
% interpMapMovieHbR = zeros([size(interpMapHbR) nIter]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure for live display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hLive = figure;
set(hLive,'color','w')
set(hLive,'name','Live display')
% Pixel units
set(hLive, 'Units', 'pixels');
screenSize = get(0,'Screensize');
screenSize = [0  300 screenSize(3) - 400 screenSize(4) - 300];
% Maximize figure
set(hLive, 'OuterPosition', screenSize);
% Colorbar Limits
hbLim = 200;
haxHbO = subplot(121); himHbO = imagesc(interpMapHbO, [-hbLim hbLim]); axis image; colorbar
title('HbO_2','FontSize',14);
% Axes set to manual
set(haxHbO, 'xlimmode','manual',...
           'ylimmode','manual',...
           'zlimmode','manual',...
           'climmode','manual',...
           'alimmode','manual');
haxHbR = subplot(122); himHbR = imagesc(interpMapHbR, [-hbLim hbLim]); axis image; colorbar
title('HbR','FontSize',14);
set(haxHbR, 'xlimmode','manual',...
           'ylimmode','manual',...
           'zlimmode','manual',...
           'climmode','manual',...
           'alimmode','manual');
% Axis position in pixels
set(haxHbO,'Units','pixels')
poshaxHbO = get(haxHbO,'Position');
widthhaxHbO = poshaxHbO(3);
heighthaxHbO = poshaxHbO(4);
set(haxHbO,'Position',[poshaxHbO(1) poshaxHbO(2) size(interpMapHbO,2) size(interpMapHbO,1)]);
set(haxHbR,'Units','pixels')
poshaxHbR = get(haxHbR,'Position');
widthhaxHbR = poshaxHbR(3);
heighthaxHbR = poshaxHbR(4);
set(haxHbR,'Position',[poshaxHbR(1) poshaxHbR(2) size(interpMapHbR,2) size(interpMapHbR,1)]);

% linkdata on
drawnow();
set(gcf,'doublebuffer','off');

%% Kalman estimation
tic
for k = 2:nIter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Estimated signal
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Time update
    [xhatminus(:, k), Pminus(:, k)] = nirs_kalman_time_update...
        (xhat(:, k-1), P(:, k-1), Q);
    % Measurement update
    [K(:, k), xhat(:, k), P(:, k) ] = nirs_kalman_measurement_update...
        (Pminus(:, k) , R, xhatminus(:, k), z(:, k));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Estimation of DC value
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Time update
    [xhatminus_DC(:, k), Pminus_DC(:, k)] = nirs_kalman_time_update...
        (xhat_DC(:, k-1), P_DC(:, k-1), Q_DC);
    % Measurement update
    [K_DC(:, k), xhat_DC(:, k), P_DC(:, k) ] = nirs_kalman_measurement_update...
        (Pminus_DC(:, k) , R_DC, xhatminus_DC(:, k), z(:, k));
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Optical Density computation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    OD(:, k) = real(log10(xhat(:, k) ./ xhat_DC(:, k)));
    % Multiply by 1e6 to get micromolar units negative sign so that an increase
    % in chromophore concentration corresponds to a decrease in intensity due to
    % light absorption
    OD(:, k) = -1e6 * OD(:, k) ./ EPF2(k, :)'; %PP used to be d = d ./ EPF2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Conversion to HbO & HbR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Modified Beer-Lambert Law
    % MBLL - c consists now of HbO and HbR, even if we had more than two
    % wavelengths to begin
    c(:, k) = inv_exs2 * OD(:, k);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Topographical projection
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % hbO = c(:, k);
    % Choose HbO channels, wl = 1
    hbO = c(NIRS.Cf.H.C.wl == 1, k);
    % Choose only channels that are visible from this view
    hbO = hbO(W.ch);
    % hbR = c(:, k);
    % Choose HbR channels, wl = 2
    hbR = c(NIRS.Cf.H.C.wl == 2, k);
    % Choose only channels that are visible from this view
    hbR = hbR(W.ch);
    % Call nirs_interpolation_render_compute as many times as needed in the loop
    interpMapHbO = nirs_interpolation_render_compute(Qinterp, hbO, interpMap);
    interpMapHbR = nirs_interpolation_render_compute(Qinterp, hbR, interpMap);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Live display
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(himHbO, 'CData', interpMapHbO);
    set(himHbR, 'CData', interpMapHbR);
    drawnow();
end
eTime = toc;
fprintf('Average processing time = %.4g\n', eTime/(nIter - 1)); 
% ~13ms no display, only computation
% ~25-32ms with live display

%% Plot figures
h2 = figure;
set(h2, 'color', 'w')
set(h2, 'Name', 'Kalman estimations')
subplot(211);
hold on
plot(t, z, 'r+', 'MarkerSize', 3);                  % noisy measurements
plot(t, xhat, 'b-', 'LineWidth',2);                 % a posteriori estimate
plot(t, xhat_DC, 'm--', 'LineWidth',2);             % a posteriori estimate (DC)
plot(eeg_t, eegAmp*sz_vector, 'k--', 'LineWidth',2);% seizure onset

xlim([t(2) t(end)])
title(sprintf('R = %0.4g, R_{DC} = %0.4g', R(1,1), R_DC(1,1)), 'FontSize', 12)
legend({'noisy measurements' 'a posteriori estimate' 'DC a posteriori estimate' 'seizure onset'},...
    'Location', 'SouthEast')
set(gca , 'FontSize', 12)
xlabel('t [s]', 'FontSize', 12); ylabel('Output', 'FontSize', 12); 

% Plot optical density
subplot(212);
hold on
plot(t(2:end), OD(:, 2:end), 'k-', 'LineWidth',2);
xlim([t(2) t(end)])
set(gca , 'FontSize', 12)
xlabel('t [s]', 'FontSize', 12); ylabel('O.D.', 'FontSize', 12); 

% Plot error covariance
figure(h1); subplot(212);
plot(t(2:end), Pminus(:, 2:end), 'k-', 'LineWidth',2);
title(sprintf('R = %0.4g', R(1,1)), 'FontSize', 12)
legend('a priori error estimate')
set(gca , 'FontSize', 12)
xlabel('t [s]', 'FontSize', 12); ylabel('(Output)^2', 'FontSize', 12); 
xlim([t(1) t(end)])

% Plot Hb concentrations
h3 = figure;
set(h3, 'color', 'w')
set(h3, 'Name', 'Hemoglobin concentrations')
subplot(211);
hold on
plot(t, c(NIRS.Cf.H.C.wl == 1,:), 'r-', 'LineWidth',1);             % HbO
plot(t, c(NIRS.Cf.H.C.wl == 2,:), 'b-', 'LineWidth',1);             % HbR
plot(eeg_t, eegAmp*sz_vector, 'k--', 'LineWidth',2);                  % seizure onset
xlim([t(2) t(end)])
title(sprintf('Kalman-estimated Hemoglobin concentrations'), 'FontSize', 12)
legend({'HbO_2' 'HbR' },...
    'Location', 'SouthEast')
set(gca , 'FontSize', 12)
xlabel('t [s]', 'FontSize', 12); ylabel('Hb [\muM]', 'FontSize', 12);

% Read offline processed data 
dataNIRS = fopen_NIR(NIRS.Dt.fir.pp(4).p{iFiles}, NIRS.Cf.H.C.N);
subplot(212)
hold on
plot(t, dataNIRS(NIRS.Cf.H.C.wl == 1,:), 'r-', 'LineWidth',1);      % HbO
plot(t, dataNIRS(NIRS.Cf.H.C.wl == 2,:), 'b-', 'LineWidth',1);      % HbR
plot(eeg_t, eegAmp*sz_vector, 'k--', 'LineWidth',2);                  % seizure onset
xlim([t(2) t(end)])
title(sprintf('Offline processing'), 'FontSize', 12)
legend({'HbO_2' 'HbR' },...
    'Location', 'SouthEast')
set(gca , 'FontSize', 12)
xlabel('t [s]', 'FontSize', 12); ylabel('Hb [\muM]', 'FontSize', 12); 

% EOF
