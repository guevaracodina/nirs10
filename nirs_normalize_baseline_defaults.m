function nirs_normalize_baseline_defaults
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________
global nirs10

nirs10.preprocessNIRS.normalize_baseline.Normalize_OD = 0; %median: 0
nirs10.preprocessNIRS.normalize_baseline.Analyzer_sf = 1; %100
nirs10.preprocessNIRS.normalize_baseline.add_or_mult = 0;
nirs10.preprocessNIRS.normalize_baseline.normalization_type = 2; %by intervals free of bad points
nirs10.preprocessNIRS.normalize_baseline.baseline_duration = 2; %2 seconds