function nirs_criugm_paces_defaults
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________
global nirs10
nirs10.preprocessNIRS.criugm_paces1.remove_no_heartbeat = 1;
%nirs10.preprocessNIRS.criugm_paces1.heart_rate_cfg = {heart_resting};
%configuration for resting state
nirs10.preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_resting.STFT_param.win_type = 0;
nirs10.preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_resting.STFT_param.win_width = 20;
nirs10.preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_resting.STFT_param.Nprobe = 200;
nirs10.preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_resting.STFT_param.fft_size = 1024;
%nirs10.preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_resting.detect_wavelength = 1;
nirs10.preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_resting.MinHeartRate = 45;
nirs10.preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_resting.MaxHeartRate = 120;
nirs10.preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_resting.InternalMinHeartRate = 45;
nirs10.preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_resting.InternalMaxHeartRate = 120; %or 3?
nirs10.preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_resting.MaxHeartStdev = 10;
%configuration for aerobic exercise
nirs10.preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_exercise.STFT_param2.win_type2 = 0;
nirs10.preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_exercise.STFT_param2.win_width2 = 20;
nirs10.preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_exercise.STFT_param2.Nprobe2 = 200;
nirs10.preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_exercise.STFT_param2.fft_size2 = 1024;
% nirs10.preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_exercise.MinHeartRate2 = 0.5;
% nirs10.preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_exercise.MaxHeartRate2 = 2;
% nirs10.preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_exercise.InternalMinHeartRate2 = 0.5;
% nirs10.preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_exercise.InternalMaxHeartRate2 = 5; %or 3?
% nirs10.preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_exercise.MaxHeartStdev2 = 0.3;
