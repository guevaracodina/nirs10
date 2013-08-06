function prior_choice = nirs_dfg_prior_choice
prior_choice        = cfg_menu;
prior_choice.name   = 'Choose priors';
prior_choice.tag    = 'prior_choice';
prior_choice.labels = {'default priors','finger_tapping_priors'};
prior_choice.values = {0,1};
prior_choice.val    = {0};
prior_choice.help   = {'Choose priors'}';