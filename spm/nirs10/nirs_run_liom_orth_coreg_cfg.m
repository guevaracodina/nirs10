function liom_OrthCoreg = nirs_run_liom_orth_coreg_cfg
%Configuration file for nirs_run_liom_orth_coreg function
%Ke Peng, 
%2013-07-25, version 0.1, 

[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'OrthCoreg');

render_image         = cfg_files;
render_image.tag     = 'render_image';
render_image.name    = 'Render image (Optional)';
render_image.help    = {'Select an image for rendering on. the T1 image will be rendered by default.'};
render_image.ufilter = '.*';
render_image.val{1}  = {''};
render_image.num     = [0 1];

channel_optode_select = cfg_menu;
channel_optode_select.tag  = 'channel_optode_select';
channel_optode_select.name = 'Select channel or optode type';
channel_optode_select.labels = {'NIRS channel' 'NIRS source' 'NIRS detector'};
channel_optode_select.values = {0 1 2};
channel_optode_select.val  = {0};
channel_optode_select.help = {'Choose results for which hemoglobin.'
                       'Must select from these three options:'
                       'HbT/HbR/HbO'}';      
                   
hemoglobin_select = cfg_menu;
hemoglobin_select.tag  = 'hemoglobin_select';
hemoglobin_select.name = 'Select hemoglobin type';
hemoglobin_select.labels = {'HbT' 'HbR' 'HbO'};
hemoglobin_select.values = {0 1 2};
hemoglobin_select.val  = {1};
hemoglobin_select.help = {'Choose results for which hemoglobin.'
                       'Must select from these three options:'
                       'HbT/HbR/HbO'}';                                   

session_select      = cfg_entry;
session_select.tag  = 'contrasts_session';
session_select.name = 'Input session number';
session_select.val = {0};
session_select.strtype = 'r';
session_select.num     = [1 1];
session_select.help    = {'Specify sessions. Input the number. Input 0 for a group view result.'};  
                   
correction_select = cfg_menu;
correction_select.tag  = 'correction_select';
correction_select.name = 'Select correction type';
correction_select.labels = {'Uncorrected' 'EC corrected'};
correction_select.values = {0 1};
correction_select.val  = {0};
correction_select.help = {'Choose results with which correction method.'
                       'Must select from these two options:'
                       'Uncorrected/EC corrected'}';                   
                   
NIRS_channels_optodes         = cfg_branch;
NIRS_channels_optodes.tag     = 'NIRS_channels_optodes';
NIRS_channels_optodes.name    = 'Coregister NIRS channels or optodes';
NIRS_channels_optodes.val     = {channel_optode_select};
NIRS_channels_optodes.help    = {'Coregister NIRS channel positions or optode positions'};

NIRS_activation         = cfg_branch;
NIRS_activation.tag     = 'NIRS_activation';
NIRS_activation.name    = 'Coregister activation results';
NIRS_activation.val     = {hemoglobin_select session_select correction_select};
NIRS_activation.help    = {'Coregister activations'};
           
Coreg_type        = cfg_choice;
Coreg_type.name   = 'Choose Coregistration type';
Coreg_type.tag    = 'coreg_type';
Coreg_type.values = {NIRS_channels_optodes, NIRS_activation};
Coreg_type.val    = {NIRS_channels_optodes};
Coreg_type.help   = {'Choose whether to coregister NIRS optodes or to coregister activation.'}';

Coreg_layer = cfg_menu;
Coreg_layer.tag  = 'coreg_layer';
Coreg_layer.name = 'Select coregistration layer';
Coreg_layer.labels = {'scalp' 'cortex'};
Coreg_layer.values = {0 1};
Coreg_layer.val  = {0};
Coreg_layer.help = {'Choose results on which layer.'
                       'Must select from these two options:'
                       'scalp/cortex'}';   

%Excecutable branch
                   
liom_OrthCoreg      = cfg_exbranch;
liom_OrthCoreg.name = 'LIOM Orthogonal Coregistration';
liom_OrthCoreg.tag  = 'liom_OrthCoreg';
liom_OrthCoreg.val  = {NIRSmat redo1 NIRSmatCopyChoice render_image Coreg_type Coreg_layer};
liom_OrthCoreg.prog = @nirs_run_liom_orth_coreg;
liom_OrthCoreg.vout = @nirs_cfg_vout_liom_orth_coreg;
liom_OrthCoreg.help = {'This module can now be run by itself or as part of a larger batch.'}';

function vout = nirs_cfg_vout_liom_orth_coreg(job)
vout = cfg_dep;
vout.sname      = 'SPM.mat';
vout.src_output = substruct('.','SPMmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
