function nirs_normalize_baseline_defaults
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________
global nirs10

nirs10.preprocessNIRS.normalize_baseline.Normalize_OD = 0; %median: 0
nirs10.preprocessNIRS.normalize_baseline.Analyzer_sf = 1; %100
nirs10.preprocessNIRS.normalize_baseline.add_or_mult = 1;