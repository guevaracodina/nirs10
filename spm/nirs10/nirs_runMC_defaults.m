function nirs_runMC_defaults
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________
% Sets the defaults input files for runMC

global nirs10

nirs10.runMC1.mcx_exe = 'D:\Users\Philippe Pouliot\spm_MC_NIRS\Fred_MCX\mcx.exe';
nirs10.runMC1.tMCimg_exe = 'D:\Users\Philippe Pouliot\spm_MC_NIRS\tMCimg.exe';
%nirs10.runMC1.MC_runCUDAchoice = '0'; %'MCX1';
nirs10.runMC1.total_threads = '3584';
nirs10.runMC1.threads_per_block = '128';
nirs10.runMC1.gates = '10';
nirs10.runMC1.mphotons = '100000';
nirs10.runMC1.repetition = '10';
nirs10.runMC1.output_log = '1'; %1 to output log to a file; 0 to output to the screen
nirs10.runMC1.a_array = '0'; %Matlab array: 0; C array: 1
nirs10.runMC1.b_reflect = '0'; %1: to reflect photons at boundary; 0: to exit

