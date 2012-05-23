function nirs10 = tbx_cfg_nirs10
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________

addpath(fileparts(which(mfilename)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration main modules
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%module 1
readNIRS        = cfg_choice;
readNIRS.name   = 'Read NIRS data';
readNIRS.tag    = 'readNIRS';
readNIRS.values = {nirs_run_boxy_cfg nirs_run_criugm_cfg nirs_run_boxy_manual_cfg};
readNIRS.help   = {'These modules read NIRS data in different formats.'};

%module 0: utilities to read onsets and create GLM stimuli structure
readOnsets        = cfg_choice;
readOnsets.name   = 'Read Onsets';
readOnsets.tag    = 'readOnsets';
readOnsets.values = {nirs_run_AnalyzerOnsets_cfg nirs_run_readEprimeOnsets_cfg ...
    nirs_run_permuteOnsets_cfg nirs_run_addTestStimuli_cfg};
readOnsets.help   = {'These modules create stimuli structures '
    'as inputs to the General Linear Model.'}';

%module 2
preprocANAT        = cfg_choice;
preprocANAT.name   = 'Preprocess anatomical image';
preprocANAT.tag    = 'preprocANAT';
preprocANAT.values = {nirs_run_detectVitamins_cfg nirs_run_MCsegment_cfg...
    nirs_run_buildroi_cfg};
preprocANAT.help   = {'These modules pre-process anatomical images '
    'so that clean anatomical images can be used with functional data.'}';

%module 3
coregNIRS        = cfg_choice; %cfg_repeat;
coregNIRS.name   = 'Coregister NIRS data';
coregNIRS.tag    = 'coregNIRS';
coregNIRS.values = {nirs_run_coreg_cfg nirs_run_coreg_2templateT1_cfg ...
    nirs_run_coreg_helmtemp_cfg nirs_run_coreg_manual_cfg nirs_run_view3d_cfg};
coregNIRS.help   = {'These modules perform coregistration ',...
    'between NIRS and an anatomical image.'};

%module 4 - NIRS preprocessing (heart rate detection, pruning bad channels, filters)
preprocessNIRS        = cfg_choice;
preprocessNIRS.name   = 'Preprocess NIRS data';
preprocessNIRS.tag    = 'preprocessNIRS';
preprocessNIRS.values = {nirs_run_remove_chn_stdev_cfg ...
    nirs_run_criugm_paces_cfg nirs_run_mark_movement_cfg ...
     nirs_run_normalize_baseline_cfg nirs_run_ODtoHbOHbR_cfg ...
     nirs_run_generate_vhdr_vmrk_cfg}; %mark_negative HPF_LPF
preprocessNIRS.help   = {'These modules preprocess NIRS data: '
    'heart rate check, '
    'downsampling, removal of bad channels, filters.'}';

%module 8 - reconstruction
model_reconstruct    = cfg_choice; %cfg_repeat;
model_reconstruct.name   = '3D Reconstruction of NIRS data';
model_reconstruct.tag    = 'model_reconstruct';
model_reconstruct.values = {nirs_run_inverse_tikhonov_cfg nirs_run_ReMLreconstruct_cfg ...
    nirs_run_checkreconstruct_cfg nirs_run_testreconstruct_cfg};
model_reconstruct.help   = {'3D Reconstruction of NIRS data.'};

%module 9
model_specify        = cfg_choice; %cfg_repeat;
model_specify.name   = 'GLM Specification (or Averaging)';
model_specify.tag    = 'model_specify';
model_specify.values = {nirs_run_liom_intrasubject_average_cfg ...
                        nirs_run_liom_GLM_specify_cfg};
model_specify.help   = {'These modules either perform averaging or specify a GLM.'};

%module 10
model_estimate        = cfg_choice; %cfg_repeat;
model_estimate.name   = 'GLM Estimation';
model_estimate.tag    = 'model_estimate';
model_estimate.values = {nirs_run_liom_GLM_estimate_cfg nirs_run_liom_contrast_cfg  ...
    nirs_run_liom_group_cfg nirs_run_extract_map_data_cfg nirs_run_liom_1way_anova_cfg ...
    nirs_run_liom_2way_anova_cfg nirs_run_AnalyzeGLM_cfg nirs_run_ROCtest_cfg}; %liom_2way_anova
model_estimate.help   = {'These modules estimate a GLM.'};

%module 13
CRIUGM        = cfg_choice; %cfg_repeat;
CRIUGM.name   = 'CRIUGM';
CRIUGM.tag    = 'CRIUGM';
CRIUGM.values = {nirs_run_runVOIRE_cfg nirs_run_runMOB_cfg nirs_run_runOvertraining_cfg};
CRIUGM.help   = {'Data analysis for CRIUGM projects'};

%module 10
nirs_utilities        = cfg_choice;
nirs_utilities.name   = 'NIRS utilities';
nirs_utilities.tag    = 'nirs_utilities';
nirs_utilities.values = {nirs_run_nii_to_2D_cfg}; 
nirs_utilities.help   = {'Various utilities.'};

%module 11
aMRI_utilities        = cfg_choice;
aMRI_utilities.name   = 'Anatomical MRI utilities';
aMRI_utilities.tag    = 'aMRI_utilities';
aMRI_utilities.values = {nirs_resize_cfg nirs_run_thickness_cfg}; 
aMRI_utilities.help   = {'Various utilities for T1 images.'};

%-----------------------------------------------------------------------
nirs10        = cfg_choice;
nirs10.name   = 'nirs10';
nirs10.tag    = 'nirs10'; %Careful, this tag nirs10 must be the same as
%the name of the toolbox and when called by spm_jobman in nirs10.m
nirs10.values = {readNIRS readOnsets preprocessNIRS preprocANAT coregNIRS aMRI_utilities...
    nirs_run_configMC2_cfg nirs_run_runMC_cfg nirs_run_generate_sensitivity_matrix_cfg...
    nirs_run_calculatePVE_cfg model_reconstruct model_specify ...
    model_estimate nirs_utilities nirs_run_NIRS_HDM_cfg nirs_run_liom_HDM_cfg CRIUGM};