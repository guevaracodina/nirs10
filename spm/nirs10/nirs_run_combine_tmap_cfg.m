function liom_ctmap = nirs_run_combine_tmap_cfg
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'Group');
display_options = liom_contrast_group_options;
consess = nirs_spm_get_consess;

OtherNIRSmat         = cfg_files; %Select NIRS.mat for this subject
OtherNIRSmat.name    = 'Other NIRS.mat to be combined'; % The displayed name
OtherNIRSmat.help    = {'Select other NIRS.mat for the subject(s).'}; % help text displayed
OtherNIRSmat.tag     = 'OtherNIRSmat';       %file names
OtherNIRSmat.filter  = 'mat';
OtherNIRSmat.ufilter = '^NIRS.mat$';
OtherNIRSmat.num     = [1 Inf];     % Number of inputs required

FFX_or_RFX = cfg_menu;
FFX_or_RFX.tag  = 'FFX_or_RFX';
FFX_or_RFX.name = 'Fixed (patient-level) or random effects(group level)';
FFX_or_RFX.labels = {'FFX','RFX'};
FFX_or_RFX.values = {1,0};
FFX_or_RFX.val = {1};
FFX_or_RFX.help = {'Use fixed effects (FFX) for group of sessions (intra-subject) '
    'Use random effects (RFX) for group of subjects (inter-subject)'
    'FFX amounts to setting the between session variance to zero.'
    'Several subjects can be specified for FFX; they will each be treated separately.'
    'RFX assumes there is only one session per subject.'
    'RFX can also be used for one subject with multiple sessions, '
    'to take into account the variance between sessions.'}';

StatMethod      = cfg_menu;
StatMethod.tag  = 'StatMethod';
StatMethod.name = 'Statistical method for spatial correlations';
StatMethod.labels = {'EC/LKC','Bonf','2D peak-FDR','None'};
StatMethod.values = {1,2,3,0};
StatMethod.val  = {1};
StatMethod.help = {'Choose statistical method to account '
    'for false positives due to spatial correlations.'
    '(Family-wise error rate)'
    'Preferred choice: Euler characteristic calculated via Lipschitz-Killing curvature'
    'Other choice: Bonferroni correction/2D peak-FDR correction.'}';

%%%%%%%%%%%%%%%%%%%%%%%%%
group_session_to_average         = cfg_entry;
group_session_to_average.name    = 'Session to average';
group_session_to_average.tag     = 'group_session_to_average';
group_session_to_average.strtype = 'r';
group_session_to_average.num     = [1 1];
group_session_to_average.val     = {1};
group_session_to_average.help    = {'This is only used for multi-session group studies.'
    'Specify here which session the group analysis should be done upon. '
    'Only one session can be specified.'}';

group_dir_name        = cfg_entry;
group_dir_name.name    = 'Name of folder to store the analysis';
group_dir_name.tag     = 'group_dir_name';
group_dir_name.strtype = 's';
group_dir_name.num     = [1 Inf];
group_dir_name.val{1}  = 'Group';
group_dir_name.help    = {'Enter name of folder to store the analysis.'}';


number_dir_to_remove         = cfg_entry;
number_dir_to_remove.name    = 'Number of directories to remove to specify Group data storage';
number_dir_to_remove.tag     = 'number_dir_to_remove';
number_dir_to_remove.strtype = 'r';
number_dir_to_remove.num     = [1 1];
number_dir_to_remove.val     = {3};
number_dir_to_remove.help    = {'Specify the location where group data will be saved.'
    'This will be calculated by removing the specified number of directories'
    'from the file path of the first subject NIRS.mat file.'
    'Typically, 3 folders need to be removed, but sometimes only 1 or 2, sometimes 4 or more.'}';

simple_sum = cfg_menu;
simple_sum.tag    = 'simple_sum';
simple_sum.name   = 'Use simple average or precision-weighted average';
simple_sum.labels = {'Simple average','Precision-weighted average'};
simple_sum.values = {1,0};
simple_sum.val = {0};
simple_sum.help   = {'Choice of averaging at the group level.'
    'Simple average is the most conservative'
    'Precision-weighting may give more activations, but they are'
    'more likely to be false positives.'}';

contrast_p_value         = cfg_entry;
contrast_p_value.name    = 'Corrected p_value';
contrast_p_value.tag     = 'contrast_p_value';
contrast_p_value.strtype = 'r';
contrast_p_value.num     = [1 1];
contrast_p_value.val     = {0.05};
contrast_p_value.help    = {'Corrected p_value'};

user_contrasts         = cfg_branch;
user_contrasts.tag     = 'user_contrasts';
user_contrasts.name    = 'User-defined contrasts'; 
user_contrasts.val     = {consess};
user_contrasts.help    = {'User-defined contrasts.'}';

automated_contrasts         = cfg_branch;
automated_contrasts.tag     = 'automated_contrasts';
automated_contrasts.name    = 'Automated contrasts'; 
automated_contrasts.val     = {};
automated_contrasts.help    = {'Automated contrasts.'}';

ContrastChoice           = cfg_choice;
ContrastChoice.name      = 'Choose contrast choice method';
ContrastChoice.tag       = 'ContrastChoice';
ContrastChoice.values    = {user_contrasts automated_contrasts}; 
ContrastChoice.val       = {automated_contrasts}; 
ContrastChoice.help      = {'Choose method to generate contrasts'}'; 

two_sample_t_test = nirs_dfg_two_sample_t_test;
ignore_brain_mask = nirs_dfg_ignore_brain_mask;

save_big_topo = cfg_menu;
save_big_topo.tag  = 'save_big_topo';
save_big_topo.name = 'save_big_TOPO';
save_big_topo.labels = {'Yes','No'};
save_big_topo.values = {1,0};
save_big_topo.val = {1};
save_big_topo.help = {'ignore brain mask. This is used to potentially '
    'increase the size of the overlap common to all subjects.'}';

% Executable Branch
liom_ctmap      = cfg_exbranch;
liom_ctmap.name = 'Liom Combine Tmap';
liom_ctmap.tag  = 'liom_ctmap';
liom_ctmap.val  = {NIRSmat redo1 NIRSmatCopyChoice ...
    group_dir_name number_dir_to_remove OtherNIRSmat FFX_or_RFX ContrastChoice ...
    StatMethod contrast_p_value two_sample_t_test ignore_brain_mask save_big_topo ...
     group_session_to_average simple_sum display_options}; % factorial_design};
liom_ctmap.prog = @nirs_run_combine_tmap;
liom_ctmap.vout = @nirs_cfg_vout_ctmap;
liom_ctmap.help = {'Liom Group level model estimation.'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vout = nirs_cfg_vout_ctmap(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
