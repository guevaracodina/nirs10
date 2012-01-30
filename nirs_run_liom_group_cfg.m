function liom_group = nirs_run_liom_group_cfg
[NIRSmat NIRSmatCopyChoice] = get_common_NIRSmat(1,'Group');
display_options = liom_contrast_group_options;
consess = spm_get_consess;
%factorial_design = spm_get_factorial_design;

% session_number         = cfg_entry;
% session_number.name    = 'Session number';
% session_number.tag     = 'session_number';
% session_number.strtype = 'r';
% session_number.num     = [1 1];
% session_number.help    = {'Enter the number of the session you want to analyse.'
%     'Only one session can be analysed at a time.'}';

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
StatMethod.labels = {'EC/LKC','Bonf'};
StatMethod.values = {1,0};
StatMethod.val  = {1};
StatMethod.help = {'Choose statistical method to account '
    'for false positives due to spatial correlations.'
    '(Family-wise error rate)'
    'Preferred choice: Euler characteristic calculated via Lipschitz-Killing curvature'
    'Other choice: Bonferroni correction.'}';


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
contrast_p_value.name    = 'Contrast uncorrected p_value';
contrast_p_value.tag     = 'contrast_p_value';
contrast_p_value.strtype = 'r';
contrast_p_value.num     = [1 1];
contrast_p_value.val     = {0.05};
contrast_p_value.help    = {'Contrast uncorrected p_value'};

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

% Executable Branch
liom_group      = cfg_exbranch;
liom_group.name = 'Liom Group Model Estimation';
liom_group.tag  = 'liom_group';
liom_group.val  = {NIRSmat ...
    group_dir_name FFX_or_RFX ContrastChoice ...
    StatMethod contrast_p_value ...
     group_session_to_average simple_sum display_options}; % factorial_design};
liom_group.prog = @nirs_run_liom_group;
liom_group.vout = @nirs_cfg_vout_liom_group;
liom_group.help = {'Liom Group level model estimation.'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vout = nirs_cfg_vout_liom_group(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
