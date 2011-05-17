function nirs_mark_movement_defaults
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________
global nirs10

nirs10.preprocessNIRS.mark_movement.mvt_window_length = 1; %in seconds 
nirs10.preprocessNIRS.mark_movement.mvt_cutoff = 25; %as a percentage of
%signal intensity over the window length
nirs10.preprocessNIRS.mark_movement.sum_mvt_threshold = 20; %written as a
%percentage, i.e. = 10 means 10% or 0.1.
nirs10.preprocessNIRS.mark_movement.min_session_duration = 60; %in seconds
nirs10.preprocessNIRS.mark_movement.mvt_ch_thresh = 50; %percentage