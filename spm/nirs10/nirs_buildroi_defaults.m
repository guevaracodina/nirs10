function nirs_buildroi_defaults
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________
% Sets the defaults for NIRS

global nirs10

% Build ROI Options
%=======================================================================

%- output default prefix
nirs10.preprocANAT.buildroi1.output_prefix = 'roi_';
% nirs10.preprocANAT.buildroi1.lby = 'Coreg';