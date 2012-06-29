function hdm_display_options = nirs_dfg_hdm_display_options
only_display      = cfg_menu;
only_display.tag  = 'only_display';
only_display.name = 'Only display';
only_display.labels = {'Yes','No'};
only_display.values = {1,0};
only_display.val  = {0};
only_display.help = {'Only display:'
    'Use this option if this HDM has already been generated.'
    'With this option, it will not be regenerated, but display options can be changed.'
    'It will overwrite figures, so if you want to keep the old figures, you should rename'
    'the old figure directory, but not the location where the HDM.mat structure is located.'}';

show_normalized_parameters      = cfg_menu;
show_normalized_parameters.tag  = 'show_normalized_parameters';
show_normalized_parameters.name = 'show normalized parameters';
show_normalized_parameters.labels = {'Yes','No'};
show_normalized_parameters.values = {1,0};
show_normalized_parameters.val  = {0};
show_normalized_parameters.help = {'Show normalized values parameters.'}';

show_mse      = cfg_menu;
show_mse.tag  = 'show_mse';
show_mse.name = 'Show MSE';
show_mse.labels = {'Yes','No'};
show_mse.values = {1,0};
show_mse.val  = {1};
show_mse.help = {'Show mean square error on figures.'}';

plot_algebraic_CMRO2      = cfg_menu;
plot_algebraic_CMRO2.tag  = 'plot_algebraic_CMRO2';
plot_algebraic_CMRO2.name = 'Plot algebraic CMRO2';
plot_algebraic_CMRO2.labels = {'Yes','No'};
plot_algebraic_CMRO2.values = {1,0};
plot_algebraic_CMRO2.val  = {1};
plot_algebraic_CMRO2.help = {'Plot algebraic CMRO2.'}';

hdm_display_options         = cfg_branch;
hdm_display_options.tag     = 'hdm_display_options';
hdm_display_options.name    = 'HDM display options';
hdm_display_options.val     = {only_display show_normalized_parameters show_mse plot_algebraic_CMRO2}; 
hdm_display_options.help    = {'Choose HDM display options.'};
