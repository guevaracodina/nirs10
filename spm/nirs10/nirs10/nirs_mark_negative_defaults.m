function nirs_mark_negative_defaults
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________

global nirs10

nirs10.preprocessNIRS.mark_negative.sum_neg_threshold = 10; %written as a
%percentage, i.e. = 10 means 10% or 0.1.
