function ignore_brain_mask = nirs_dfg_ignore_brain_mask
ignore_brain_mask = cfg_menu;
ignore_brain_mask.tag  = 'ignore_brain_mask';
ignore_brain_mask.name = 'Ignore brain mask';
ignore_brain_mask.labels = {'Yes','No'};
ignore_brain_mask.values = {1,0};
ignore_brain_mask.val = {0};
ignore_brain_mask.help = {'ignore brain mask. This is used to potentially '
    'increase the size of the overlap common to all subjects.'}';
