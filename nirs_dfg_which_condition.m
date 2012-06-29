function which_condition = nirs_dfg_which_condition
which_condition     = cfg_entry;
which_condition.name    = 'Which condition(s)?';
which_condition.tag     = 'which_condition';
which_condition.strtype = 'r';
which_condition.val     = {1};
which_condition.num     = [0 Inf];
which_condition.help    = {'Enter stimuli numbers to include as a Matlab row '
    'vector (get from the design matrix associated with the data file)'}';
