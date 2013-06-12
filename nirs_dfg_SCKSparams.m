function SCKSparams = nirs_dfg_SCKSparams
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCKS parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SCKSnoise      = cfg_menu;
SCKSnoise.tag  = 'SCKSnoise';
SCKSnoise.name = 'Promote parameters and noise to states';
SCKSnoise.labels = {'Yes','No'};
SCKSnoise.values = {1,0};
SCKSnoise.val  = {0};
SCKSnoise.help = {'Promote parameters and noise to time-dependent states.'}';
    
State_annealing         = cfg_entry; 
State_annealing.name    = 'Annealing factor on states';
State_annealing.tag     = 'State_annealing';       
State_annealing.strtype = 'r';
State_annealing.num     = [1 1];     
State_annealing.val     = {0.9995};
State_annealing.help    = {'Enter annealing factor on states. This is a number between 0 and 1.'
    'A value closer to 1 corresponds to more annealing. A value of 0.5 corresponds to very little annealing.'}';

Parameter_annealing         = cfg_entry; 
Parameter_annealing.name    = 'Annealing factor on parameters';
Parameter_annealing.tag     = 'Parameter_annealing';       
Parameter_annealing.strtype = 'r';
Parameter_annealing.num     = [1 1];     
Parameter_annealing.val     = {0.9995};
Parameter_annealing.help    = {'Enter annealing factor on parameters.' 
    'This is a number between 0 and 1.'
    'A value closer to 1 corresponds to more annealing. A value of 0.5 corresponds to very little annealing.'}';

SCKSparams         = cfg_branch;
SCKSparams.tag     = 'SCKSparams';
SCKSparams.name    = 'SCKS parameters';
SCKSparams.val     = {SCKSnoise State_annealing Parameter_annealing}; 
SCKSparams.help    = {'User-controlled SCKS parameters.'};

