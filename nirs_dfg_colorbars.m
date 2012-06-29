function [override_colorbar GroupColorbars] = nirs_dfg_colorbars
GroupColorbars      = cfg_menu;
GroupColorbars.tag  = 'GroupColorbars';
GroupColorbars.name = 'Group or separate colorbars';
GroupColorbars.labels = {'Group','Separate'};
GroupColorbars.values = {1,0};
GroupColorbars.val = {0};
GroupColorbars.help = {'When considering deactivations (inverted responses),'
    'This allows displaying one (group) or two (separate) colorbars'}';


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
