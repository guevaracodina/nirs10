function nirs_MCsegment_defaults
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________
% Sets the defaults for NIRS

global nirs10

% Monte Carlo Segmentation Options
%=======================================================================

% MC segmentation parameters
nirs10.preprocANAT.MCsegment1.output_autonaming = 0;
nirs10.preprocANAT.MCsegment1.output_prefix = 'Only edit if you chose ''No'' to ''Automatic output naming'''; % output default name
nirs10.preprocANAT.MCsegment1.thresh_as = 0.6; % sorting threshold for voxel with no conflict in belonging 
nirs10.preprocANAT.MCsegment1.rebel_surrounding = 3; % surrounding size for rebel voxels
nirs10.preprocANAT.MCsegment1.rebel_thresh_hs = 0.3; % threshold for head shadow
%- sorting method
nirs10.preprocANAT.MCsegment1.wtm = 0; % White Matter
nirs10.preprocANAT.MCsegment1.grm = 0; % Grey Matter
nirs10.preprocANAT.MCsegment1.csf = 0; % CSF
nirs10.preprocANAT.MCsegment1.skl = 2; % Skull
nirs10.preprocANAT.MCsegment1.skn = 1; % Skin

% head shadow parameters
nirs10.preprocANAT.head_shadow.se_size_hs = 2; % size of the structuring element (Mathematical morphology)
nirs10.preprocANAT.head_shadow.thresh_hs = 0.6; % threshold in mc_segment_hs to get a uniformised mask of the head

% process image parameters
nirs10.preprocANAT.process_image.se_size_pi = 2;
nirs10.preprocANAT.process_image.gaussfilt_size = 7;
nirs10.preprocANAT.process_image.gaussfilt_sdev = 4;