function config = nirs_interpolation_render_config (brain_view, AllowExtrapolation, no_interpolation, new_path, figure_name, thz)
% SYNTAX
% config = nirs_interpolation_render_config (brain_view, AllowExtrapolation,
% no_interpolation, new_path, figure_name, thz)
% INPUTS
% brain_view            % Brain view
%                       1 % 'ventral'
%                       2 % 'dorsal'
%                       3 %'right_lateral'
%                       4 %'left_lateral'
%                       5 %'frontal'
%                       6 %'occipital'
% AllowExtrapolation    0: do not extrapolate
% no_interpolation      0: interpolate
% new_path              Path where interpolated images are saved
% figure_name           Name of interpolated image
% thz                   Threshold value of the interpolated image. 
%                       Attention: Cannot be zero
% OUTPUT
% config                Configuration structure
%_______________________________________________________________________________
% Copyright (C) 2014 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

% Brain view, or a loop can be used here to generate all views
config.brain_view = brain_view; 
% Display chosen brain view 
[~, spec_hemi] = nirs_get_brain_view(config.brain_view);
fprintf('Brain view: %s\n', spec_hemi);
% Extrapolation
config.AllowExtrapolation = AllowExtrapolation; 
% Interpolation
config.no_interpolation = no_interpolation; 
% Path to save image
config.new_path = new_path; 
% Figure name
config.figure_name = figure_name; 
% Threshold value
config.thz = thz; 

% EOF
