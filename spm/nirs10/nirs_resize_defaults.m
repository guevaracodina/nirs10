function nirs_resize_defaults
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________
% Sets the defaults input files for runMC

global nirs10

nirs10.coregNIRS.resize1.out_autonaming =0;
nirs10.coregNIRS.resize1.out_dim =[1 1 1];
nirs10.coregNIRS.resize1.out_dt ='same';