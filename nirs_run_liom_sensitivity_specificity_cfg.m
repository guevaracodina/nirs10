function liom_sensitivity_specificity = nirs_run_liom_sensitivity_specificity_cfg
%Configuration file for nirs_run_liom_projection function
%Ke Peng, 
%2012-08-08, version 0.1, 
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'Sen_spe');

focus_location      = cfg_entry;
focus_location.tag  = 'focus_location';
focus_location.name = 'focus location (MNI)';
focus_location.val{1} = [19 52 26];
focus_location.strtype = 'r';
focus_location.num     = [3 Inf];
focus_location.help    = {'Enter focus_location in mm in MNI coordinates.'
    'Make sure this is consistent with the anatomical image for this subject.'
    'Specify a list of focus locations, with one row for the 3 coordinates of each subject'}';

focus_radius      = cfg_entry;
focus_radius.tag  = 'focus_radius';
focus_radius.name = 'focus radius in mm';
focus_radius.val = {6};
focus_radius.strtype = 'r';
focus_radius.num     = [1 Inf];
focus_radius.help    = {'input focus radius'
    'If the focus radii are different between subjects, '
    'Specify a list of focus radii, one for each subject'}';

%Style

render_style      = cfg_menu;
render_style.tag  = 'render_style';
render_style.name = 'Style';
render_style.labels = {'new','old'};
render_style.values = {0,1};
render_style.val  = {0};
render_style.help = {'Choose render style.'
                       'Must select from these two options:'
                       'new/old'}';

%Brighten blobs

render_blobs      = cfg_menu;
render_blobs.tag  = 'render_blobs';
render_blobs.name = 'Brighten blobs';
render_blobs.labels = {'none','slightly','more','lots'};
render_blobs.values = {0,1,2,3};
render_blobs.val  = {0};
render_blobs.help = {'Choose brighten blobs style.'
                       'Must select from these four options:'
                       'none/slightly/more/lots'}';

%which colours

render_colour      = cfg_menu;
render_colour.tag  = 'render_colour';
render_colour.name = 'Which colours';
render_colour.labels = {'RGB'};% "custom" option has been removed
render_colour.values = {0};
render_colour.val  = {0};
render_colour.help = {'Choose display colours.'
                       'Must select from these two options:'
                       'RGB/Custom'}';
                   
%Sensitivity and specificity config
sessions_sel      = cfg_menu;
sessions_sel.tag  = 'sessions_sel';
sessions_sel.name = 'Sessions to check for each subject';
sessions_sel.labels = {'All sessions'};
sessions_sel.values = {0};
sessions_sel.val  = {0};
sessions_sel.help = {'For one single subject, select sessions that we look at to check sensitivity/specificity'
                       'Must select from these options:'
                       'All sessions'}';

thres_sel      = cfg_menu;
thres_sel.tag  = 'thres_sel';
thres_sel.name = 'Apply threshold value for constrasts';
thres_sel.labels = {'No','Yes'};
thres_sel.values = {0,1};
thres_sel.val  = {1};
thres_sel.help = {'Use the uncorrected t-map or use the corrected t-map'
                       'Must select from these options:'
                       'No/Yes'}';
                   
%Excecutable branch
                   
liom_sensitivity_specificity      = cfg_exbranch;
liom_sensitivity_specificity.name = 'LIOM Sensitivity Specificity Check';
liom_sensitivity_specificity.tag  = 'liom_sensitivity_specificity';
liom_sensitivity_specificity.val  = {NIRSmat redo1 NIRSmatCopyChoice focus_location ...
    focus_radius render_style render_blobs render_colour sessions_sel thres_sel};
liom_sensitivity_specificity.prog = @nirs_run_liom_sensitivity_specificity;
liom_sensitivity_specificity.vout = @nirs_cfg_vout_liom_sensitivity_specificity;
liom_sensitivity_specificity.help = {'This module can now be run by itself or as part of a larger batch.'}';

function vout = nirs_cfg_vout_liom_sensitivity_specificity(job)
vout = cfg_dep;
vout.sname      = 'SPM.mat';
vout.src_output = substruct('.','SPMmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});


