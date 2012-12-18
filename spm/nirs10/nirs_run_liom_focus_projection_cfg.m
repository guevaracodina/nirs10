function liom_projection = nirs_run_liom_focus_projection_cfg
%Configuration file for nirs_run_liom_projection function
%Ke Peng, 
%2012-08-08, version 0.1, 
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'Proj');

focus_location      = cfg_entry;
focus_location.tag  = 'focus_location';
focus_location.name = 'focus location (MNI)';
focus_location.val = {[19 52 26]};
focus_location.strtype = 'r';
focus_location.num     = [Inf 3];
focus_location.help    = {'Enter focus_location in mm in MNI coordinates.'
    'Make sure this is consistent with the anotomical image for this subject.'
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

%Project contrasts configurations                   
                   
contrasts_views      = cfg_entry;
contrasts_views.tag  = 'contrasts_views';
contrasts_views.name = 'Select Views to generate';
contrasts_views.val = {[2 3 4 5]};
contrasts_views.strtype = 'r';
contrasts_views.num     = [1 inf];
contrasts_views.help    = {'Specify views. Input the corresponding number.'
                           '1-Ventral;2-Dorsal;3-Right;4-Left;5-Frontal;6-Occipital'}';                     

contrasts_session      = cfg_entry;
contrasts_session.tag  = 'contrasts_session';
contrasts_session.name = 'Select Sessions to project contrasts';
contrasts_session.val = {1};
contrasts_session.strtype = 'r';
contrasts_session.num     = [1 inf];
contrasts_session.help    = {'Specify sessions. Input the number.'};                   
                   
contrasts_disabled         = cfg_branch;
contrasts_disabled.tag     = 'contrasts_disabled';
contrasts_disabled.name    = 'Do not project contrasts';
contrasts_disabled.val     = {};
contrasts_disabled.help    = {'Disable the contrasts projection'};
                                     
contrasts_enabled         = cfg_branch;
contrasts_enabled.tag     = 'contrasts_enabled';
contrasts_enabled.name    = 'Project contrasts onto a same image';
contrasts_enabled.val     = {contrasts_session contrasts_views};
contrasts_enabled.help    = {'Enable the contrasts projection'};

proj_contrasts      = cfg_choice;
proj_contrasts.tag  = 'proj_contrasts';
proj_contrasts.name = 'Contrasts Projection';
proj_contrasts.values = {contrasts_enabled contrasts_disabled};
proj_contrasts.val = {contrasts_disabled};
proj_contrasts.help = {'Choose whether to project contrasts on the same image.'
                       'Must select from these two options:'
                       'Projuct contrasts onto a same image/Do not project contrasts'}';

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
                   
liom_projection      = cfg_exbranch;
liom_projection.name = 'LIOM Contrasts Projection';
liom_projection.tag  = 'liom_projection';
liom_projection.val  = {NIRSmat redo1 NIRSmatCopyChoice focus_location ...
    focus_radius render_style render_blobs render_colour proj_contrasts thres_sel};
liom_projection.prog = @nirs_run_liom_focus_projection;
liom_projection.vout = @nirs_cfg_vout_liom_focus_projection;
liom_projection.help = {'This module can now be run by itself or as part of a larger batch.'}';

function vout = nirs_cfg_vout_liom_focus_projection(job)
vout = cfg_dep;
vout.sname      = 'SPM.mat';
vout.src_output = substruct('.','SPMmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});

