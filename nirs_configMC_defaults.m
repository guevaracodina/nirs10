function nirs_configMC_defaults
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________
% Sets the defaults input files for configMC

global nirs10

nirs10.configMC1.MC_CUDAchoice = 1;
nirs10.configMC1.MC_configdir = 'MC'; %no longer used

nirs10.configMC1.nphotons = 1e8;
nirs10.configMC1.seed = 394581069;
nirs10.configMC1.modulationFreq = 0;
nirs10.configMC1.deltaT = 8e-9;
nirs10.configMC1.numTimeGates = 1;
nirs10.configMC1.radiis = 0.4;
nirs10.configMC1.radiid = 2.5;
nirs10.configMC1.voxelSize = 1;
%Note that I inverted the order of wavelengths compared to the
%code of Michèle

nirs10.configMC1.perturbationPpties_l1 = [0.0000   0.00   0.0   0.0];

nirs10.configMC1.perturbationPpties_l2 = [0.0000   0.00   0.0   0.0];
