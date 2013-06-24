function display_options = liom_cine_options
contrast_figures      = cfg_menu;
contrast_figures.tag  = 'contrast_figures';
contrast_figures.name = 'Generate figures';
contrast_figures.labels = {'No','Both .fig and .png','Only .fig','Only .png'};
contrast_figures.values = {0,1,2,3};
contrast_figures.val = {0};
contrast_figures.help = {'Generate contrast figures. '
    'Note .fig colorbar is incorrect - it is not saved properly by Matlab.'
    'Use .png to view colorbar for t-stat.'}';

override_colorbar = nirs_dfg_colorbars_cine;
show_boundary = nirs_dfg_show_boundary;

display_options         = cfg_branch;
display_options.tag     = 'display_options';
display_options.name    = 'Display options'; 
display_options.val     = {contrast_figures override_colorbar show_boundary};
display_options.help    = {'Display options or rarely used options.'}';
