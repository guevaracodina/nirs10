function nirs_paces_defaults
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________
global nirs10

nirs10.preprocessNIRS.paces1.win_type = 0;
nirs10.preprocessNIRS.paces1.win_width = 20;
nirs10.preprocessNIRS.paces1.Nprobe = 200;
nirs10.preprocessNIRS.paces1.fft_size = 1024;
nirs10.preprocessNIRS.paces1.remove_no_heartbeat = 1;
nirs10.preprocessNIRS.paces1.detect_wavelength = 1;
nirs10.preprocessNIRS.paces1.MinHeartRate = 0.5;
nirs10.preprocessNIRS.paces1.MaxHeartRate = 2;
nirs10.preprocessNIRS.paces1.InternalMinHeartRate = 0.5;
nirs10.preprocessNIRS.paces1.InternalMaxHeartRate = 5; %or 3?
nirs10.preprocessNIRS.paces1.MaxHeartStdev = 0.3;