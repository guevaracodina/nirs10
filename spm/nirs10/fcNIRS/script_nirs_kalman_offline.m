% Script to read NIRS data and perform Kalman filtering offline
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
preProcStep = 4;
nSessions = numel(NIRS.Dt.fir.pp(preProcStep).p);
% Choose file(session)
iFiles = 3;
dataNIRS = fopen_NIR(NIRS.Dt.fir.pp(preProcStep).p{iFiles}, NIRS.Cf.H.C.N);
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
round(eeg_fs*(t(end) - t(1)))
eeg_t = linspace(t(1), t(end), round(eeg_fs*(t(end) - t(1))));
sz_vector = zeros(size(eeg_t));
sz_start = [107704 ];
sz_end = [160302 ];
sz_vector(sz_start(1):sz_end(1)) = 1;
figure; stem(eeg_t, sz_vector)

%% Plot data
h1 = figure; set(h1,'color','w')
set(h1,'Name',sprintf('Pre-processing step: %s', NIRS.Dt.fir.pp(preProcStep).pre))
subplot(211)
fprintf('Wavelength: %dnm & %dnm\n',NIRS.Cf.dev.wl)
hb = 1;
switch hb
    case 1
        hbLabel = 'HbO';
    case 2
        hbLabel = 'HbR';
    case 3
        hbLabel = 'HbT';
end
plot(t, dataNIRS(NIRS.Cf.H.C.wl == hb,:)');
ylabel('[\muM]')
xlabel('t [s]')
% title(['\lambda' sprintf(' = %d nm',NIRS.Cf.dev.wl(hb))])
title(hbLabel)

%% Kalman Example (Time-Varying Design)
% Consider the discrete plant with additive Gaussian noise w[n] on the input
% u[n] and data:
% 
% x[n+1] = Ax[n] + B(u[n] +w[n])
% y[n]   = Cx[n] 

A = [1.1269   -0.4940    0.1129
     1.0000         0         0
          0    1.0000         0];

B = [-0.3832
      0.5919
      0.5191];

C = [1 0 0]; 

% Our goal is to design a Kalman filter that estimates the output y[n] given the
% inputs u[n] and the noisy output measurements
% y_v[n] = Cx[n] + v[n]
% where v[n] is some Gaussian white noise. 

% Assuming that Q = R = 1, you can now design the discrete Kalman filter by

Q = 0.1; R = 100;

% Generate a sinusoidal input u and process and measurement noise vectors w and v.
% t = (0:100)';
% u = sin(t/5);
u = mean(dataNIRS(NIRS.Cf.H.C.wl == hb,:))';

n = length(t);
randn('seed',0);
w = sqrt(Q)*randn(n,1);
v = sqrt(R)*randn(n,1);

% Use process noise w and measurement noise v generated above
sys = ss(A,B,C,0,-1);
y = lsim(sys,u+w);      % w = process noise
yv = y + v;             % v = measurement noise

%% Implement the time-varying filter with the following for loop, given the
% initial conditions

P = B*Q*B';         % Initial error covariance
x = zeros(3,1);     % Initial condition on the state
ye = zeros(length(t),1);
ycov = zeros(length(t),1); 

for i=1:length(t)
  % Measurement update
  Mn = P*C'/(C*P*C'+R);
  x = x + Mn*(yv(i)-C*x);   % x[n|n]
  P = (eye(3)-Mn*C)*P;      % P[n|n]

  ye(i) = C*x;
  errcov(i) = C*P*C';

  % Time update
  x = A*x + B*u(i);        % x[n+1|n]
  P = A*P*A' + B*Q*B';     % P[n+1|n]
end

%% You can now compare the true and estimated output graphically.
h2 = figure;
set(h2,'color','w')
set(h2,'Name', 'Kalman estimation')
subplot(211), 
plot(t,y/100,'r.',t,ye/100,'k-', eeg_t, sz_vector, 'b--')
legend({'True response' 'Filtered response' 'Seizure onset'})
title(['Time-varying Kalman filter response (' hbLabel ')'])
xlabel('t [s]'), ylabel(['\Delta' hbLabel '/' hbLabel '(%)'])
subplot(212), plot(t,y-yv,'-.',t,y-ye,'-')
legend({'measurement error' 'estimation error'})
xlabel('t [s]'), ylabel('Output')
% The first plot shows the true response y (dashed line) and the filtered
% response ye (solid line). The second plot compares the measurement error
% (dash-dot) with the estimation error (solid).

%%
% The time-varying filter also estimates the covariance errcov of the estimation
% error y ? ye at each sample. Plot it to see if your filter reached steady
% state (as you expect with stationary input noise).
figure(h1);
subplot(212)
plot(t,errcov), ylabel('Error covar'), xlabel('t [s]')

EstErr = y-ye;
EstErrCov = sum(EstErr.*EstErr)/length(EstErr)

%% New figure
h2 = figure;
set(h2,'color','w')
set(h2,'Name', 'Seizure')
hb = 1; % HbO
switch hb
    case 1
        hbLabel = 'HbO';
    case 2
        hbLabel = 'HbR';
    case 3
        hbLabel = 'HbT';
end
y = mean(dataNIRS(NIRS.Cf.H.C.wl == hb,:))';
hold on
plot(t,y/100,'r-','LineWidth',2)
% e = (std(dataNIRS(NIRS.Cf.H.C.wl == hb,:))')/sqrt(size(dataNIRS,1));
% downSamplingFactor = 100;
% errorbar(t(1:downSamplingFactor:end),y(1:downSamplingFactor:end)/100,e(1:downSamplingFactor:end),'r-')
plot(eeg_t, sz_vector, 'k--','LineWidth',2)
hb = 2; % HbR
switch hb
    case 1
        hbLabel = 'HbO';
    case 2
        hbLabel = 'HbR';
    case 3
        hbLabel = 'HbT';
end
y = mean(dataNIRS(NIRS.Cf.H.C.wl == hb,:))';
plot(t,y/100,'b-','LineWidth',2)
set(gca,'FontSize',14)
legend({'HbO' 'Seizure onset' 'HbR'})
title('Average response','FontSize',14)
xlabel('t [s]','FontSize',14), ylabel('\DeltaHb/Hb [%]','FontSize',14)

% EOF
