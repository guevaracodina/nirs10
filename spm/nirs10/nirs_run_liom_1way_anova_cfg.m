function liom_1way_anova = nirs_run_liom_1way_anova_cfg
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'OneWayAnova');
display_options = liom_contrast_group_options;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% One-way anova
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
anova_level         = cfg_entry;
anova_level.name    = 'Number of levels';
anova_level.tag     = 'anova_level';
anova_level.strtype = 'r';
anova_level.num     = [1 1];
anova_level.val     = {2};
anova_level.help    = {'Enter number of levels (e.g. 2 for young vs old).'
    'For a one-way anova, there is only one factor.'}';

level_name        = cfg_entry;
level_name.name    = 'Name for this level';
level_name.tag     = 'level_name';
level_name.strtype = 's';
level_name.num     = [1 Inf];
level_name.help    = {'Enter name of this level.'}';

level_subj        = cfg_entry;
level_subj.name    = 'Subjects at this level';
level_subj.tag     = 'level_subj';
level_subj.strtype = 'r';
level_subj.num     = [1 Inf];
level_subj.help    = {'Enter subject numbers for this level.'}';

level         = cfg_branch;
level.tag     = 'level';
level.name    = 'Add level';
level.val     = {level_name level_subj};
level.help    = {'Add level'}';

level_repeat         = cfg_repeat;
level_repeat.tag     = 'level_repeat';
level_repeat.name    = 'New level';
level_repeat.help    = {'Add a level for this factor'}';
level_repeat.values  = {level};
level_repeat.num     = [1 Inf];

anova_dir_name        = cfg_entry;
anova_dir_name.name    = 'Name of folder to store the analysis';
anova_dir_name.tag     = 'anova_dir_name';
anova_dir_name.strtype = 's';
anova_dir_name.num     = [1 Inf];
anova_dir_name.val{1}  = 'Anova';
anova_dir_name.help    = {'Enter name of folder to store the analysis.'}';

contrast_p_value         = cfg_entry;
contrast_p_value.name    = 'Corrected p_value';
contrast_p_value.tag     = 'contrast_p_value';
contrast_p_value.strtype = 'r';
contrast_p_value.num     = [1 1];
contrast_p_value.val     = {0.05};
contrast_p_value.help    = {'Corrected p_value'};

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

% Executable Branch
liom_1way_anova      = cfg_exbranch;
liom_1way_anova.name = 'Liom 1-way Anova Estimation';
liom_1way_anova.tag  = 'liom_1way_anova';
liom_1way_anova.val  = {NIRSmat redo1 NIRSmatCopyChoice ...
    anova_dir_name number_dir_to_remove anova_level level_repeat ...
    StatMethod contrast_p_value ...
    group_session_to_average display_options}; % factorial_design};
liom_1way_anova.prog = @nirs_run_liom_1way_anova;
liom_1way_anova.vout = @nirs_cfg_vout_liom_1way_anova;
liom_1way_anova.help = {'Liom 1way anova estimation.'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vout = nirs_cfg_vout_liom_1way_anova(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});