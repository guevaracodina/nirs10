function EM_parameters = nirs_dfg_hdm_EM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EM parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Niterations         = cfg_entry; 
Niterations.name    = 'Maximum number of EM iterations';
Niterations.tag     = 'Niterations';       
Niterations.strtype = 'r';
Niterations.num     = [1 1];     
Niterations.val     = {128};
Niterations.help    = {'Maximum number of EM iterations. 128 is the basic number.'
    'Increase to 512 or more if necessary.'}';

kernel_window         = cfg_entry; 
kernel_window.name    = 'Kernel window length in seconds';
kernel_window.tag     = 'kernel_window';       
kernel_window.strtype = 'r';
kernel_window.num     = [1 1];     
kernel_window.val     = {20};
kernel_window.help    = {'Kernel window length in seconds.'}';

dFcriterion         = cfg_entry; 
dFcriterion.name    = 'Convergence criterion';
dFcriterion.tag     = 'dFcriterion';       
dFcriterion.strtype = 'r';
dFcriterion.num     = [1 1];     
dFcriterion.val     = {1};
dFcriterion.help    = {'Convergence criterion on changes of free energy F.'
    'Changes in F less than this value are required for convergence in E and in M steps.'
    '1e-2 is the SPM8 default value.'}';

LogAscentRate         = cfg_entry; 
LogAscentRate.name    = 'Initial log ascent rate';
LogAscentRate.tag     = 'LogAscentRate';       
LogAscentRate.strtype = 'r';
LogAscentRate.num     = [1 1];     
LogAscentRate.val     = {-2};
LogAscentRate.help    = {'Initial log ascent rate: control initial rate of movement in parameter space'}';

MaxLogAscentRate         = cfg_entry; 
MaxLogAscentRate.name    = 'Maximal log ascent rate';
MaxLogAscentRate.tag     = 'MaxLogAscentRate';       
MaxLogAscentRate.strtype = 'r';
MaxLogAscentRate.num     = [1 1];     
MaxLogAscentRate.val     = {4};
MaxLogAscentRate.help    = {'Maximal absolute value of log ascent rate: control minimal/maximal rate of movement in parameter space'}';

spm_integrator      = cfg_menu;
spm_integrator.tag  = 'spm_integrator';
spm_integrator.name = 'Choose ODE integrator';
spm_integrator.labels = {'spm_int','spm_int_ode','spm_int_J'};
spm_integrator.values = {'spm_int','spm_int_ode','spm_int_J'};
spm_integrator.val  = {'spm_int'};
spm_integrator.help = {'Choose integrator to use for ordinary differential equations.'
    'spm_int is fastest, spm_int_ode most precise, spm_int_J in between'}';

Mstep_iterations         = cfg_entry; 
Mstep_iterations.name    = 'Maximum number of M step iterations';
Mstep_iterations.tag     = 'Mstep_iterations';       
Mstep_iterations.strtype = 'r';
Mstep_iterations.num     = [1 1];     
Mstep_iterations.val     = {8};
Mstep_iterations.help    = {'Maximum number of M-step iterations. 8 is the standard choice.'}';

EM_parameters         = cfg_branch;
EM_parameters.tag     = 'EM_parameters';
EM_parameters.name    = 'Parameters of EM algorithm';
EM_parameters.val     = {Niterations kernel_window spm_integrator ... 
    dFcriterion LogAscentRate MaxLogAscentRate Mstep_iterations}; 
EM_parameters.help    = {'Parameters of Expectation Maximization algorithm.'};
