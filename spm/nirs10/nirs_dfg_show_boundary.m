function show_boundary = nirs_dfg_show_boundary
show_boundary      = cfg_menu;
show_boundary.tag  = 'show_boundary';
show_boundary.name = 'Show mask boundary';
show_boundary.labels = {'Yes','No'};
show_boundary.values = {1,0};
show_boundary.val  = {1};
show_boundary.help = {'This allows displaying a fine black line along the boundary of the mask.'}';
