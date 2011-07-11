function nirs_criugm_defaults
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________
% Sets the defaults for NIRS coregistration

global nirs10

nirs10.readNIRS.criugm1.generic2.subj.age1 = 25;
nirs10.readNIRS.criugm1.heart_pace = 1;
nirs10.readNIRS.criugm1.modify_prop_helmet = 0;
nirs10.readNIRS.criugm1.baseline_method = 2;
nirs10.readNIRS.criugm1.CWsystem = 6;

nirs10.readNIRS.criugm1.generic2.subj.bs_pos.nasion = 'Nasion';
nirs10.readNIRS.criugm1.generic2.subj.bs_pos.leftear = 'LeftEar';
nirs10.readNIRS.criugm1.generic2.subj.bs_pos.rightear = 'RightEar';

