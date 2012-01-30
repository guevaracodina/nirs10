function display_options = liom_contrast_group_options

contrast_figures      = cfg_menu;
contrast_figures.tag  = 'contrast_figures';
contrast_figures.name = 'Generate figures';
contrast_figures.labels = {'No','Both .fig and .tiff','Only .fig','Only .tiff'};
contrast_figures.values = {0,1,2,3};
contrast_figures.val = {3};
contrast_figures.help = {'Generate contrast figures. '
    'Note .fig colorbar is incorrect - it is not saved properly by Matlab.'
    'Use .tiff to view colorbar for t-stat.'}';

figures_visible      = cfg_menu;
figures_visible.tag  = 'figures_visible';
figures_visible.name = 'Make figures visible';
figures_visible.labels = {'Yes','No'};
figures_visible.values = {1,0};
figures_visible.val = {0};
figures_visible.help = {'Make figures visible during processing, and staying open at the end.'}';

colorbar_max         = cfg_entry;
colorbar_max.name    = 'Colorbar maximum value';
colorbar_max.tag     = 'colorbar_max';
colorbar_max.strtype = 'r';
colorbar_max.num     = [1 1];
colorbar_max.val     = {5};
colorbar_max.help    = {'Enter maximum value for colorbar'};

colorbar_min         = cfg_entry;
colorbar_min.name    = 'Colorbar minimum value';
colorbar_min.tag     = 'colorbar_min';
colorbar_min.strtype = 'r';
colorbar_min.num     = [1 1];
colorbar_min.val     = {2};
colorbar_min.help    = {'Enter minimum value for colorbar'};

colorbar_max2         = cfg_entry;
colorbar_max2.name    = 'Colorbar maximum value';
colorbar_max2.tag     = 'colorbar_max2';
colorbar_max2.strtype = 'r';
colorbar_max2.num     = [1 1];
colorbar_max2.val     = {-2};
colorbar_max2.help    = {'Enter maximum value for colorbar for negative maps'};

colorbar_min2         = cfg_entry;
colorbar_min2.name    = 'Colorbar minimum value';
colorbar_min2.tag     = 'colorbar_min2';
colorbar_min2.strtype = 'r';
colorbar_min2.num     = [1 1];
colorbar_min2.val     = {-5};
colorbar_min2.help    = {'Enter minimum value for colorbar for negative maps'};

colorbar_override      = cfg_branch;
colorbar_override.name      = 'Override colorbar';
colorbar_override.tag       = 'colorbar_override';
colorbar_override.val       = {colorbar_min colorbar_max colorbar_min2 colorbar_max2};
colorbar_override.help      = {'Override colorbar.'};

colorbar_default      = cfg_branch;
colorbar_default.name      = 'Default colorbar';
colorbar_default.tag       = 'colorbar_default';
colorbar_default.val       = {};
colorbar_default.help      = {'Default colorbar.'};

override_colorbar           = cfg_choice;
override_colorbar.name      = 'Override colorbar';
override_colorbar.tag       = 'override_colorbar';
override_colorbar.values    = {colorbar_default colorbar_override};
override_colorbar.val       = {colorbar_default};
override_colorbar.help      = {'Override default treatment of colorbar.'
    'User can then specify maximum and minimum values for the colorbar.'}';

GenerateInverted      = cfg_menu;
GenerateInverted.tag  = 'GenerateInverted';
GenerateInverted.name = 'Generate Inverted Responses';
GenerateInverted.labels = {'Yes','No'};
GenerateInverted.values = {1,0};
GenerateInverted.val = {1};
GenerateInverted.help = {'Generate contrasts for inverted responses.'};

GroupColorbars      = cfg_menu;
GroupColorbars.tag  = 'GroupColorbars';
GroupColorbars.name = 'Group or separate colorbars';
GroupColorbars.labels = {'Group','Separate'};
GroupColorbars.values = {1,0};
GroupColorbars.val = {0};
GroupColorbars.help = {'When considering deactivations (inverted responses),'
    'This allows displaying one (group) or two (separate) colorbars'}';

SmallFigures      = cfg_menu;
SmallFigures.tag  = 'SmallFigures';
SmallFigures.name = 'Large or small figures';
SmallFigures.labels = {'Large','Small'};
SmallFigures.values = {0,1};
SmallFigures.val = {1};
SmallFigures.help = {'Write to disk large or small (compressed) figures.'}';

GroupFiguresIntoSubplots      = cfg_menu;
GroupFiguresIntoSubplots.tag  = 'GroupFiguresIntoSubplots';
GroupFiguresIntoSubplots.name = 'Group Figures Into Subplots';
GroupFiguresIntoSubplots.labels = {'Yes','No'};
GroupFiguresIntoSubplots.values = {1,0};
GroupFiguresIntoSubplots.val = {1};
GroupFiguresIntoSubplots.help = {'Group Figures Into Subplots.'};

output_unc      = cfg_menu;
output_unc.tag  = 'output_unc';
output_unc.name = 'Output unc. figures';
output_unc.labels = {'Yes','No'};
output_unc.values = {1,0};
output_unc.val = {0};
output_unc.help = {'Output figures that are uncorrected against false positives.'}';

write_neg_pos      = cfg_menu;
write_neg_pos.tag  = 'write_neg_pos';
write_neg_pos.name = 'Write neg/pos separate contrasts';
write_neg_pos.labels = {'Yes', 'No'};
write_neg_pos.values = {1,0};
write_neg_pos.val    = {0};
write_neg_pos.help = {'If generating negative contrasts, whether to output '
    'separate maps for negative and positive contrasts and for both, '
    'or only the maps with both contrasts' }';

save_nifti_contrasts      = cfg_menu;
save_nifti_contrasts.tag  = 'save_nifti_contrasts';
save_nifti_contrasts.name = 'Save contrast images in nifti format';
save_nifti_contrasts.labels = {'Yes', 'No'};
save_nifti_contrasts.values = {1,0};
save_nifti_contrasts.val    = {0};
save_nifti_contrasts.help = {'This option is useful for 2nd level studies: '
    'It creates nifti images, one for each contrast, which can then be input'
    'Into an SPM 2nd level analysis.' }';

display_options         = cfg_branch;
display_options.tag     = 'display_options';
display_options.name    = 'Display options'; 
display_options.val     = {GenerateInverted GroupColorbars ...
    contrast_figures override_colorbar figures_visible GroupFiguresIntoSubplots ...
    output_unc SmallFigures write_neg_pos save_nifti_contrasts};
display_options.help    = {'Display options or rarely used options.'}';
