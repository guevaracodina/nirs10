function nirs_boxy_defaults
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________
% Sets the defaults for NIRS

global nirs10
nirs10.readNIRS.boxy1.generic1.subj.age1 = 25;
nirs10.readNIRS.boxy1.config_path.prj_path = 'mtg';
nirs10.readNIRS.boxy1.config_path.T1_path = 'T1';
nirs10.readNIRS.boxy1.config_path.output_path = 'dataSPM';
nirs10.readNIRS.boxy1.cf1.save_bin1 = true;
nirs10.readNIRS.boxy1.cf1.Lambda = [830 690]; %Wavelengths (must match the indexing order)
%nirs10.readNIRS.boxy1.cf1.LambdaHbO = [1 0]; 
nirs10.readNIRS.boxy1.cf1.freq = 19.5312;
%maximum distance between source and detector to keep channel, in centimeters
nirs10.readNIRS.boxy1.cf1.distmax = 6; 
%minimum distance
nirs10.readNIRS.boxy1.cf1.distmin = 1;
% %Size of chunks of data to be read in the main conversion loop
nirs10.readNIRS.boxy1.cf1.sizebloc = 1024; %suggested to take a power of 2 due to filtering 
nirs10.readNIRS.boxy1.cf1.nb_Mux = 32; %number of MUX (NOT the number of sources)
nirs10.readNIRS.boxy1.cf1.MaxSources = 64; %maximum number of sources
nirs10.readNIRS.boxy1.cf1.nb_Det = 16; %number of detectors
nirs10.readNIRS.boxy1.cf1.use10_10system = true; %true to use 10_10system, otherwise 10_20system
nirs10.readNIRS.boxy1.cf1.MaxElectrodes = 19; %maximum number of electrodes
nirs10.readNIRS.boxy1.cf1.resample = 1; %desampling factor; use 1 to keep all data

