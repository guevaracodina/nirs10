function liom_2way_anova = nirs_run_liom_2way_anova_cfg
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'TwoWayAnova');
display_options = liom_contrast_group_options;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two-way anova
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
anova_dir_name.val{1}  = 'TwoWayAnova';
anova_dir_name.help    = {'Enter name of folder to store the analysis.'}';

anova2_sessions         = cfg_entry;
anova2_sessions.name    = 'Sessions to include';
anova2_sessions.tag     = 'anova2_sessions';
anova2_sessions.strtype = 'r';
anova2_sessions.num     = [1 Inf];
anova2_sessions.val     = {[1 2]};
anova2_sessions.help    = {'Specify which sessions to include in the anova'}';

anova2_contrasts         = cfg_entry;
anova2_contrasts.name    = 'Contrasts to include';
anova2_contrasts.tag     = 'anova2_contrasts';
anova2_contrasts.strtype = 'r';
anova2_contrasts.num     = [1 Inf];
anova2_contrasts.val     = {[1 2]};
anova2_contrasts.help    = {'Specify which contrasts to include in the anova'}';

includeSubjectEffects      = cfg_menu;
includeSubjectEffects.tag  = 'includeSubjectEffects';
includeSubjectEffects.name = 'Include Subject Effects';
includeSubjectEffects.labels = {'Yes','No'};
includeSubjectEffects.values = {1,0};
includeSubjectEffects.val = {1};
includeSubjectEffects.help = {'Include regressors to account for intrasubject effects.'}';

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

% Executable Branch
liom_2way_anova      = cfg_exbranch;
liom_2way_anova.name = 'Liom 2-way Anova Estimation';
liom_2way_anova.tag  = 'liom_2way_anova';
liom_2way_anova.val  = {NIRSmat redo1 NIRSmatCopyChoice  ...
    anova_dir_name number_dir_to_remove StatMethod anova2_sessions anova2_contrasts ...
    includeSubjectEffects contrast_p_value display_options};
liom_2way_anova.prog = @nirs_run_liom_2way_anova;
liom_2way_anova.vout = @nirs_cfg_vout_liom_2way_anova;
liom_2way_anova.help = {'Liom 2way anova estimation.'
    'This module currently assumes that the first factor comes from the list of sessions'
    'while the second factor comes from the list of contrasts.'
    'It does not allow distinguishing between different types of subjects.'}';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vout = nirs_cfg_vout_liom_2way_anova(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
