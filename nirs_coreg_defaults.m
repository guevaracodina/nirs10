function nirs_coreg_defaults
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________
% Sets the defaults for NIRS coregistration

global nirs10

nirs10.coregNIRS.coreg1.anatT1_template = {fullfile(fileparts(which('spm')),'templates','T1.nii')};
nirs10.coregNIRS.coreg1.nasion_wMNI = [0 84 -48];
nirs10.coregNIRS.coreg1.AL_wMNI = [-83 -19 -38];
nirs10.coregNIRS.coreg1.AR_wMNI = [ 83 -19 -38];
nirs10.coregNIRS.coreg1.GenDataTopo = 1;
nirs10.coregNIRS.coreg1.fid_in_subject_MNI = 0;