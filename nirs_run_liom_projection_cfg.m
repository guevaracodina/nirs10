function liom_projection = nirs_run_liom_projection_cfg
%Configuration file for nirs_run_liom_projection function
%Ke Peng, 
%2012-08-08, version 0.1, Function created
%   Detailed explanation goes here

%Select SPM

SPMmat         = cfg_files; %Select NIRS.mat for this subject
SPMmat.name    = 'SPM.mat (optional)'; % The displayed name    
SPMmat.help    = {'Select SPM.mat for the subject(s).'}; % help text displayed
SPMmat.tag     = 'SPMmat';       %file names
SPMmat.filter  = 'mat';
SPMmat.ufilter = '^SPM.mat$';
SPMmat.num     = [0 Inf];     % Number of inputs required

%Select contrasts to run

contrasts_selected      = cfg_entry;
contrasts_selected.tag  = 'contrasts_selected';
contrasts_selected.name = 'Select contrasts to run';
contrasts_selected.val = {-1};
contrasts_selected.strtype = 'r';
contrasts_selected.num     = [0 Inf];
contrasts_selected.help    = {'Specify contrasts to run. Multiple integer values could be inputted. Input -1 for all sessions and all onsets'};

%Applying mask options

mask_image_select      = cfg_files;
mask_image_select.tag  = 'mask_image_select';
mask_image_select.name = 'Select images for masking';
mask_image_select.num     = [1 1];
mask_image_select.filter  = 'nii';
mask_image_select.ufilter = {};
mask_image_select.help    = {'One image file (.nii) must be selected'};

mask_image_nature_mask      = cfg_menu;
mask_image_nature_mask.tag  = 'mask_image_nature_mask';
mask_image_nature_mask.name = 'Nature of mask';
mask_image_nature_mask.labels = {'inclusive','exclusive'};
mask_image_nature_mask.values = {0,1};
mask_image_nature_mask.val  = {0};
mask_image_nature_mask.help = {'Choose nature of mask.'
                       'Must select from these two options:'
                       'inclusive/exclusive'}';

mask_image         = cfg_branch;
mask_image.tag     = 'mask_image';
mask_image.name    = 'Image';
mask_image.val     = {mask_image_select mask_image_nature_mask};
mask_image.help    = {'Please Select image for masking'};

mask_contrast_select      = cfg_entry;
mask_contrast_select.tag  = 'mask_contrast_select';
mask_contrast_select.name = 'Select contrasts for masking';
mask_contrast_select.val = {-1};
mask_contrast_select.strtype = 'r';
mask_contrast_select.num     = [1 1];
mask_contrast_select.help    = {'Specify contrasts for masking. Input the number for certain session'};

mask_contrast_unc_p_value      = cfg_entry;
mask_contrast_unc_p_value.tag  = 'mask_contrast_unc_p_value';
mask_contrast_unc_p_value.name = 'Uncorrected mask p value';
mask_contrast_unc_p_value.val = {0.05};
mask_contrast_unc_p_value.strtype = 'r';
mask_contrast_unc_p_value.num     = [1 1];
mask_contrast_unc_p_value.help    = {'Specify the uncorrected mask p value. Real number must be inputted.'};


mask_contrast_nature_mask      = cfg_menu;
mask_contrast_nature_mask.tag  = 'mask_contrast_nature_mask';
mask_contrast_nature_mask.name = 'Nature of mask';
mask_contrast_nature_mask.labels = {'inclusive','exclusive'};
mask_contrast_nature_mask.values = {0,1};
mask_contrast_nature_mask.val  = {0};
mask_contrast_nature_mask.help = {'Choose nature of mask.'
                       'Must select from these two options:'
                       'inclusive/exclusive'}';

mask_contrast         = cfg_branch;
mask_contrast.tag     = 'mask_contrast';
mask_contrast.name    = 'Contrast';
mask_contrast.val     = {mask_contrast_select mask_contrast_unc_p_value mask_contrast_nature_mask};
mask_contrast.help    = {'Please Select contrast for masking'};

mask_none         = cfg_branch;
mask_none.tag     = 'mask_none';
mask_none.name    = 'None';
mask_none.val     = {};
mask_none.help    = {'Do not apply masks'};

mask_option      = cfg_choice;
mask_option.tag  = 'mask_option';
mask_option.name = 'Applying mask';
mask_option.values = {mask_none mask_contrast mask_image};
mask_option.val = {mask_none};
mask_option.help = {'Choose whether to apply masks.'
                       'Must select from these three options:'
                       'none/contrast/image'}';

%Title for comparison

title_comparison_input      = cfg_entry;
title_comparison_input.tag  = 'title_comparison_input';
title_comparison_input.name = 'Input comparison title';
title_comparison_input.val = {};
title_comparison_input.strtype = 'r';
title_comparison_input.num     = [1 1];
title_comparison_input.help    = {'Input your desired title for comparison.'};

title_input         = cfg_branch;
title_input.tag     = 'title_input';
title_input.name    = 'By input';
title_input.val     = {title_comparison_input};
title_input.help    = {'Title will be inputted by the operator.'};

title_default         = cfg_branch;
title_default.tag     = 'title_default';
title_default.name    = 'Default';
title_default.val     = {};
title_default.help    = {'Title will be chosen automatically according to selected contrasts.'};

title_comparison      = cfg_choice;
title_comparison.tag  = 'title_comparison';
title_comparison.name = 'Title for Comparison';
title_comparison.values = {title_default title_input};
title_comparison.val = {title_default};
title_comparison.help = {'Choose the title for comparison.'
                       'Must select from these two options:'
                       'Default/By input'}';

%p value adjusted to control

p_value_none_threshold_input      = cfg_entry;
p_value_none_threshold_input.tag  = 'p_value_none_threshold_input';
p_value_none_threshold_input.name = 'threshold (T or P value)';
p_value_none_threshold_input.val = {0.001};
p_value_none_threshold_input.strtype = 'r';
p_value_none_threshold_input.num     = [1 1];
p_value_none_threshold_input.help    = {'input threshold value'};

p_value_none         = cfg_branch;
p_value_none.tag     = 'p_value_none';
p_value_none.name    = 'None';
p_value_none.val     = {p_value_none_threshold_input};
p_value_none.help    = {'Not to use the FWE method'};

p_value_FWE_p_input      = cfg_entry;
p_value_FWE_p_input.tag  = 'p_value_FWE_p_input';
p_value_FWE_p_input.name = 'p value (FWE)';
p_value_FWE_p_input.val = {0.05};
p_value_FWE_p_input.strtype = 'r';
p_value_FWE_p_input.num     = [1 1];
p_value_FWE_p_input.help    = {'input p value for FWE method'};

p_value_FWE         = cfg_branch;
p_value_FWE.tag     = 'p_value_FWE';
p_value_FWE.name    = 'FWE';
p_value_FWE.val     = {p_value_FWE_p_input};
p_value_FWE.help    = {'To use the FWE method'};

p_value_adjustment      = cfg_choice;
p_value_adjustment.tag  = 'p_value_adjustment';
p_value_adjustment.name = 'p value adjustment to control';
p_value_adjustment.values = {p_value_FWE p_value_none};
p_value_adjustment.val = {p_value_FWE};
p_value_adjustment.help = {'Specify the p value adjustment to control.'
                       'Must select from these two options:'
                       'FWE/None'}';

%extent threshold {voxels}

extent_threshold      = cfg_entry;
extent_threshold.tag  = 'extent_threshold';
extent_threshold.name = '&extent threshold {voxels}';
extent_threshold.val = {0};
extent_threshold.strtype = 'r';
extent_threshold.num     = [1 1];
extent_threshold.help    = {'input extent threshold value'};

%Select Render file

render_file_select      = cfg_files;
render_file_select.tag  = 'render_file_select';
render_file_select.name = 'Select render file';
render_file_select.num     = [1 1];
render_file_select.filter  = 'mat';
render_file_select.ufilter = {};
render_file_select.help    = {'One render file must be selected'};

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

%Render Section

render_file_enabled         = cfg_branch;
render_file_enabled.tag     = 'render_file_enabled';
render_file_enabled.name    = 'Render to file function enabled';
render_file_enabled.val     = {render_file_select render_style render_blobs render_colour};
render_file_enabled.help    = {'Render to file activated'};

render_file_disabled         = cfg_branch;
render_file_disabled.tag     = 'render_file_disabled';
render_file_disabled.name    = 'Render to file function disabled';
render_file_disabled.val     = {};
render_file_disabled.help    = {'Render to file disactivated'};

render_file      = cfg_choice;
render_file.tag  = 'render_file';
render_file.name = 'Render to file';
render_file.values = {render_file_enabled render_file_disabled};
render_file.val = {render_file_enabled};
render_file.help = {'Specify whether render contrasts to image file'};

%excutable branch

liom_projection      = cfg_exbranch;
liom_projection.name = 'LIOM Contrasts Projection';
liom_projection.tag  = 'liom_projection';
liom_projection.val  = {SPMmat contrasts_selected mask_option title_comparison p_value_adjustment extent_threshold render_file};
liom_projection.prog = @nirs_run_liom_projection;
liom_projection.vout = @nirs_cfg_vout_liom_projection;
liom_projection.help = {'This module can now be run by itself or as part of a larger batch.'
    'Attention: Please give a correct SPM structure.'
    'At this stage, this module is not able to analyze a SPM structure generated by Ye GLM algorithm'}';

function vout = nirs_cfg_vout_liom_projection(job)
vout = cfg_dep;
vout.sname      = 'SPM.mat';
vout.src_output = substruct('.','SPMmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});

