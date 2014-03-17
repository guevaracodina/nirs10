% Script to read NIRS data and perform Kalman filtering offline
%_______________________________________________________________________________
% Copyright (C) 2014 LIOM Laboratoire d'Imagerie Optique et Moléculaire
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
% NIRSmat = 'F:\Edgar\Data\NIRS\epiNIRS_data\epiNIRS_Processed\epiNIRS\epi127SD\dataSPMa\coreg\NIRS.mat';
% epi143MDM  Good optodes coverage
% NIRSmat = 'F:\Edgar\Data\NIRS\epiNIRS_data\epiNIRS_Processed\epiNIRS\epi143MDM\dataSPMa\coreg\NIRS.mat';
% NIRSmat = 'F:\Edgar\Data\NIRS\epiNIRS_data\ControlNIRS_Processed\ControlNIRS\Sujet001\dataSPMa\coreg\NIRS.mat';
% epi122 - Seizures in sessions 
NIRSmat = 'F:\Edgar\Data\NIRS\epiNIRS_data\epiNIRS_Processed\epiNIRS\epi122PEV\dataSPMa\coreg\NIRS.mat';
load(NIRSmat)

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
%     if preProcStep == 4
    dataNIRS = fopen_NIR(NIRS.Dt.fir.pp(preProcStep).p{iFiles}, NIRS.Cf.H.C.N);
else
    dataNIRS = fopen_NIR(NIRS.Dt.fir.pp(preProcStep).p{iFiles},NIRS.Cf.H.C.id(1,end));
end
% dataNIRS = [];
% for iFiles = 1:nSessions,
%     % currentData [nChannels x nTimePoints]
%     currentData = fopen_NIR(NIRS.Dt.fir.pp(preProcStep).p{iFiles}, NIRS.Cf.H.C.N);
%     % dataNIRS [nChannels x (nTimePoints x nSessions)]
%     dataNIRS = [dataNIRS, currentData];
% end
% Time vector
t = linspace(0, (size(dataNIRS,2)-1)/NIRS.Cf.dev.fs , size(dataNIRS,2));

%% EEG seizures
% Epi122PEV (Session 1,3,4,6)
eeg_fs = 1/2e-3;
round(eeg_fs*(t(end) - t(1)));
eeg_t = linspace(t(1), t(end), round(eeg_fs*(t(end) - t(1))));
sz_vector = zeros(size(eeg_t));
sz_start = [107704 ];
sz_end = [160302 ];
sz_vector(sz_start(1):sz_end(1)) = 1;
% figure; stem(eeg_t, sz_vector)

%% Plot data
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
% z = mean(dataNIRS(NIRS.Cf.H.C.wl == hb,:));     % 1 channel
% z = dataNIRS(NIRS.Cf.H.C.wl == hb,:);           % Channels at 1 wavelength
z = dataNIRS;                                   % All Channels
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

% Initial guess the measurement is zero
x_hat_DC(:,1) = 0;
%  Initial guess Choose P_0 to be 1
P_DC(:,1) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hemoglobin concentrations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Partial volume correction factor
PVF = [50 50];
% Differential pathlength factor
% From Duncan et al., Phys. Med. Biol. 1995 40(2):295-304.
DPF = [6.51 5.86]; % FOR 690 and 832 nm
DPF = repmat(DPF',[1 nChannels]); % nLambda x nChannels
PVF = repmat(PVF',[1 nChannels]); % nLambda x nChannels
EPF = zeros(1, nChannels); %EPF = L * PPF = L * DPF/PVF at each wavelength
% Source-detector distances
Cgp = NIRS.Cf.H.C.gp;
% Channels wavelength
Cwl = NIRS.Cf.H.C.wl;
% Device wavelength
wl = NIRS.Cf.dev.wl;
for iChannels = 1:nChannels
    EPF(1,iChannels) = Cgp(1,iChannels)*DPF(Cwl(1,iChannels),iChannels)./PVF(Cwl(1,iChannels),iChannels);
end
% exs(:,1): HbO for each wavelength; exs(:,2): HbR for each wavelength Alexis's
% choice of extinction coefficients corresponds to case 1 in Homer, which
% appears to be their preferred choice too
[exs, extcoeff_ref] = GetExtinctions(wl,1);
inv_exs = pinv(exs(:,1:2)); % inv_exs has size 2 x #wl
inv_exs2 = kron(inv_exs,eye(nChannels/size(wl,2))); % size #pairs*2 x #pairs*#wl
% (number of channels nChannels = number of pairs x number of wavelengths)

%% Kalman estimation
% R is the estimate of measurement variance, change to see effect
R = 0.01*ones([nChannels 1]);
R_DC = 10*ones([nChannels 1]);

tic
for k = 2:nIter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Estimated signal
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Time update
    xhatminus(:, k)    = xhat(:, k-1);
    Pminus(:, k)       = P(:, k-1) + Q;
    
    % Measurement update
    K(:, k)            = Pminus(:, k) ./ (Pminus(:, k) + R);
    xhat(:, k)         = xhatminus(:, k) + K(:, k) .* (z(:, k) - xhatminus(:, k));
    P(:, k)            = (1 - K(:, k)) .* Pminus(:, k);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Estimation of DC value
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Time update
    xhatminus_DC(:, k)    = xhat_DC(:, k-1);
    Pminus_DC(:, k)       = P_DC(:, k-1) + Q_DC;
    
    % Measurement update
    K_DC(:, k)            = Pminus_DC(:, k) ./ (Pminus_DC(:, k) + R_DC);
    xhat_DC(:, k)         = xhatminus_DC(:, k) + K_DC(:, k) .* (z(:, k) - xhatminus_DC(:, k));
    P_DC(:, k)            = (1 - K_DC(:, k)) .* Pminus_DC(:, k);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Optical Density computation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    OD(:, k) = real(log10(xhat(:, k) ./ xhat_DC(:, k)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Conversion to HbO & HbR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Effective path length
    EPF2 = ones(size(OD,2),1)*EPF; %PP used to be EPF2 = EPF *ones(1,size(d,2));
    % Multiply by 1e6 to get micromolar units negative sign so that an increase in
    % chromophore concentration corresponds to a decrease in intensity due to light
    % absorption
    OD(:, k) = -1e6 * OD(:, k) ./ EPF2'; %PP used to be d = d ./ EPF2;
    % Modified Beer-Lambert Law
    % MBLL - c consists now of HbO and HbR, even if we had more
    % than two wavelengths to begin
    c(:, k) = inv_exs2 * OD(:, k);
  
end
eTime = toc;
fprintf('Average processing time = %.4g\n', eTime/nIter);
h = figure;
set(h, 'color', 'w')
subplot(211);
hold on
plot(t, z, 'r+', 'MarkerSize', 3);                  % noisy measurements
plot(eeg_t, 4000*sz_vector, 'k--', 'LineWidth',2);  % seizure onset
plot(t, xhat, 'b-', 'LineWidth',2);                 % a posteriori estimate
plot(t, xhat_DC, 'm--', 'LineWidth',2);             % a posteriori estimate (DC)
   
xlim([t(2) t(end)])
title(sprintf('R = %0.4g, R_{DC} = %0.4g', R(1,1), R_DC(1,1)), 'FontSize', 12)
legend({'noisy measurements' 'seizure onset' 'a posteriori estimate' 'DC a posteriori estimate'},...
    'Location', 'SouthEast')
set(gca , 'FontSize', 12)
xlabel('t [s]', 'FontSize', 12); ylabel('Output', 'FontSize', 12); 

% Plot optical density
subplot(212);
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

% %% You can now compare the true and estimated output graphically.
% h2 = figure;
% set(h2,'color','w')
% set(h2,'Name', 'Kalman estimation')
% subplot(211), 
% plot(t,y,'r.',t,ye,'k-', eeg_t, sz_vector, 'b--')
% legend({'True response' 'Filtered response' 'Seizure onset'})
% title(['Time-varying Kalman filter response (' hbLabel ')'])
% xlabel('t [s]'), ylabel(['\Delta' hbLabel '/' hbLabel '(%)'])
% subplot(212), plot(t,y-yv,'-.',t,y-ye,'-')
% legend({'measurement error' 'estimation error'})
% xlabel('t [s]'), ylabel('Output')
% % The first plot shows the true response y (dashed line) and the filtered
% % response ye (solid line). The second plot compares the measurement error
% % (dash-dot) with the estimation error (solid).
% 
% %%
% % The time-varying filter also estimates the covariance errcov of the estimation
% % error y ? ye at each sample. Plot it to see if your filter reached steady
% % state (as you expect with stationary input noise).
% figure(h1);
% subplot(212)
% plot(t,errcov), ylabel('Error covar'), xlabel('t [s]')
% 
% EstErr = y-ye;
% EstErrCov = sum(EstErr.*EstErr)/length(EstErr)
% 
% %% New figure
% h2 = figure;
% set(h2,'color','w')
% set(h2,'Name', 'Seizure')
% hb = 1; % HbO
% switch hb
%     case 1
%         hbLabel = 'HbO';
%     case 2
%         hbLabel = 'HbR';
%     case 3
%         hbLabel = 'HbT';
% end
% y = mean(dataNIRS(NIRS.Cf.H.C.wl == hb,:))';
% hold on
% plot(t,y/100,'r-','LineWidth',2)
% % e = (std(dataNIRS(NIRS.Cf.H.C.wl == hb,:))')/sqrt(size(dataNIRS,1));
% % downSamplingFactor = 100;
% % errorbar(t(1:downSamplingFactor:end),y(1:downSamplingFactor:end)/100,e(1:downSamplingFactor:end),'r-')
% plot(eeg_t, sz_vector, 'k--','LineWidth',2)
% hb = 2; % HbR
% switch hb
%     case 1
%         hbLabel = 'HbO';
%     case 2
%         hbLabel = 'HbR';
%     case 3
%         hbLabel = 'HbT';
% end
% y = mean(dataNIRS(NIRS.Cf.H.C.wl == hb,:))';
% plot(t,y/100,'b-','LineWidth',2)
% set(gca,'FontSize',14)
% legend({'HbO' 'Seizure onset' 'HbR'})
% title('Average response','FontSize',14)
% xlabel('t [s]','FontSize',14), ylabel('\DeltaHb/Hb [%]','FontSize',14)

% EOF
