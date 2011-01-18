function nirs_remove_chn_stdev_defaults
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________
global nirs10

nirs10.preprocessNIRS.remove_chn_stdev.threshold_stdev = 0.2; %written as a
%percentage, i.e. = 0.5 means 0.5% or 0.005. Values as low as 0.1 are
%reasonable
nirs10.preprocessNIRS.remove_chn_stdev.window_stdev = 5;
