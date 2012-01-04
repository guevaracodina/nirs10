function liom_1way_anova = nirs_run_liom_1way_anova_cfg
NIRSmat         = cfg_files; %Select NIRS.mat for this subject
NIRSmat.name    = 'NIRS.mat'; % The displayed name
NIRSmat.tag     = 'NIRSmat';       %file names
NIRSmat.filter  = 'mat';
NIRSmat.ufilter = '^NIRS.mat$';
NIRSmat.num     = [1 Inf];     % Number of inputs required
NIRSmat.help    = {'Select NIRS.mat for the subject(s).'}; % help text displayed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% One-way anova
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPM factorial design configuration
%%%%%%%%%%%%%%%%

% ---------------------------------------------------------------------
% dir Directory
% ---------------------------------------------------------------------
dir         = cfg_files;
dir.tag     = 'dir';
dir.name    = 'Directory';
dir.help    = {'Select a directory where the SPM.mat file containing the specified design matrix will be written.'};
dir.filter = 'dir';
dir.val{1} = {''};
dir.ufilter = '.*';
dir.num     = [0 1];

% ---------------------------------------------------------------------
% scans Scans
% ---------------------------------------------------------------------
scans         = cfg_files;
scans.tag     = 'scans';
scans.name    = 'Scans';
scans.help    = {'Select the images.  They must all have the same image dimensions, orientation, voxel size etc.'};
scans.filter = 'image';
scans.ufilter = '.*';
scans.num     = [0 Inf];
% ---------------------------------------------------------------------
% t1 One-sample t-test
% ---------------------------------------------------------------------
t1         = cfg_branch;
t1.tag     = 't1';
t1.name    = 'One-sample t-test';
t1.val     = {scans };
t1.help    = {''};

% ---------------------------------------------------------------------
% scans1 Group 1 scans
% ---------------------------------------------------------------------
scans1         = cfg_files;
scans1.tag     = 'scans1';
scans1.name    = 'Group 1 scans';
scans1.help    = {'Select the images from sample 1.  They must all have the same image dimensions, orientation, voxel size etc.'};
scans1.filter = 'image';
scans1.ufilter = '.*';
scans1.num     = [0 Inf];
% ---------------------------------------------------------------------
% scans2 Group 2 scans
% ---------------------------------------------------------------------
scans2         = cfg_files;
scans2.tag     = 'scans2';
scans2.name    = 'Group 2 scans';
scans2.help    = {'Select the images from sample 2.  They must all have the same image dimensions, orientation, voxel size etc.'};
scans2.filter = 'image';
scans2.ufilter = '.*';
scans2.num     = [0 Inf];
% ---------------------------------------------------------------------
% dept Independence
% ---------------------------------------------------------------------
dept         = cfg_menu;
dept.tag     = 'dept';
dept.name    = 'Independence';
dept.help    = {
    'By default, the measurements are assumed to be independent between levels. '
    ''
    'If you change this option to allow for dependencies, this will violate the assumption of sphericity. It would therefore be an example of non-sphericity. One such example would be where you had repeated measurements from the same subjects - it may then be the case that, over subjects, measure 1 is correlated to measure 2. '
    ''
    'Restricted Maximum Likelihood (REML): The ensuing covariance components will be estimated using ReML in spm_spm (assuming the same for all responsive voxels) and used to adjust the statistics and degrees of freedom during inference. By default spm_spm will use weighted least squares to produce Gauss-Markov or Maximum likelihood estimators using the non-sphericity structure specified at this stage. The components will be found in SPM.xVi and enter the estimation procedure exactly as the serial correlations in fMRI models.'
    ''
    }';
dept.labels  = {
    'Yes'
    'No'
    }';
dept.values  = {0 1};
dept.val     = {0};
% ---------------------------------------------------------------------
% deptn Independence (default is 'No')
% ---------------------------------------------------------------------
deptn         = cfg_menu;
deptn.tag     = 'dept';
deptn.name    = 'Independence';
deptn.help    = {
    'By default, the measurements are assumed to be dependent between levels. '
    ''
    'If you change this option to allow for dependencies, this will violate the assumption of sphericity. It would therefore be an example of non-sphericity. One such example would be where you had repeated measurements from the same subjects - it may then be the case that, over subjects, measure 1 is correlated to measure 2. '
    ''
    'Restricted Maximum Likelihood (REML): The ensuing covariance components will be estimated using ReML in spm_spm (assuming the same for all responsive voxels) and used to adjust the statistics and degrees of freedom during inference. By default spm_spm will use weighted least squares to produce Gauss-Markov or Maximum likelihood estimators using the non-sphericity structure specified at this stage. The components will be found in SPM.xVi and enter the estimation procedure exactly as the serial correlations in fMRI models.'
    ''
    }';
deptn.labels  = {
    'Yes'
    'No'
    }';
deptn.values  = {0 1};
deptn.val     = {1};

% ---------------------------------------------------------------------
% variance Variance
% ---------------------------------------------------------------------
variance         = cfg_menu;
variance.tag     = 'variance';
variance.name    = 'Variance';
variance.help    = {
    'By default, the measurements in each level are assumed to have unequal variance. '
    ''
    'This violates the assumption of ''sphericity'' and is therefore an example of ''non-sphericity''.'
    ''
    'This can occur, for example, in a 2nd-level analysis of variance, one contrast may be scaled differently from another.  Another example would be the comparison of qualitatively different dependent variables (e.g. normals vs. patients).  Different variances (heteroscedasticy) induce different error covariance components that are estimated using restricted maximum likelihood (see below).'
    ''
    'Restricted Maximum Likelihood (REML): The ensuing covariance components will be estimated using ReML in spm_spm (assuming the same for all responsive voxels) and used to adjust the statistics and degrees of freedom during inference. By default spm_spm will use weighted least squares to produce Gauss-Markov or Maximum likelihood estimators using the non-sphericity structure specified at this stage. The components will be found in SPM.xVi and enter the estimation procedure exactly as the serial correlations in fMRI models.'
    ''
    }';
variance.labels = {
    'Equal'
    'Unequal'
    }';
variance.values = {0 1};
variance.val    = {1};
% ---------------------------------------------------------------------
% gmsca Grand mean scaling
% ---------------------------------------------------------------------
gmsca         = cfg_menu;
gmsca.tag     = 'gmsca';
gmsca.name    = 'Grand mean scaling';
gmsca.help    = {
    'This option is only used for PET data.'
    ''
    'Selecting YES will specify ''grand mean scaling by factor'' which could be eg. ''grand mean scaling by subject'' if the factor is ''subject''. '
    ''
    'Since differences between subjects may be due to gain and sensitivity effects, AnCova by subject could be combined with "grand mean scaling by subject" to obtain a combination of between subject proportional scaling and within subject AnCova. '
    ''
    }';
gmsca.labels = {
    'No'
    'Yes'
    }';
gmsca.values = {0 1};
gmsca.val    = {0};
% ---------------------------------------------------------------------
% ancova ANCOVA
% ---------------------------------------------------------------------
ancova         = cfg_menu;
ancova.tag     = 'ancova';
ancova.name    = 'ANCOVA';
ancova.help    = {
    'This option is only used for PET data.'
    ''
    'Selecting YES will specify ''ANCOVA-by-factor'' regressors. This includes eg. ''Ancova by subject'' or ''Ancova by effect''. These options allow eg. different subjects to have different relationships between local and global measurements. '
    ''
    }';
ancova.labels = {
    'No'
    'Yes'
    }';
ancova.values = {0 1};
ancova.val    = {0};
% ---------------------------------------------------------------------
% t2 Two-sample t-test
% ---------------------------------------------------------------------
t2         = cfg_branch;
t2.tag     = 't2';
t2.name    = 'Two-sample t-test';
t2.val     = {scans1 scans2 dept variance gmsca ancova };
t2.help    = {''};

% ---------------------------------------------------------------------
% scans Scans [1,2]
% ---------------------------------------------------------------------
scans         = cfg_files;
scans.tag     = 'scans';
scans.name    = 'Scans [1,2]';
scans.help    = {'Select the pair of images. '};
scans.filter = 'image';
scans.ufilter = '.*';
scans.num     = [2 2];
% ---------------------------------------------------------------------
% pair Pair
% ---------------------------------------------------------------------
pair         = cfg_branch;
pair.tag     = 'pair';
pair.name    = 'Pair';
pair.val     = {scans };
pair.help    = {'Add a new pair of scans to your experimental design'};
% ---------------------------------------------------------------------
% generic Pairs
% ---------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Pairs';
generic.help    = {''};
generic.values  = {pair};
generic.num     = [1 Inf];
% ---------------------------------------------------------------------
% pt Paired t-test
% ---------------------------------------------------------------------
pt         = cfg_branch;
pt.tag     = 'pt';
pt.name    = 'Paired t-test';
pt.val     = {generic gmsca ancova};
pt.help    = {''};

% ---------------------------------------------------------------------
% scans Scans
% ---------------------------------------------------------------------
scans         = cfg_files;
scans.tag     = 'scans';
scans.name    = 'Scans';
scans.help    = {'Select the images.  They must all have the same image dimensions, orientation, voxel size etc.'};
scans.filter = 'image';
scans.ufilter = '.*';
scans.num     = [1 Inf];
% ---------------------------------------------------------------------
% c Vector
% ---------------------------------------------------------------------
c         = cfg_entry;
c.tag     = 'c';
c.name    = 'Vector';
c.help    = {'Vector of covariate values'};
c.strtype = 'e';
c.num     = [Inf 1];
% ---------------------------------------------------------------------
% cname Name
% ---------------------------------------------------------------------
cname         = cfg_entry;
cname.tag     = 'cname';
cname.name    = 'Name';
cname.help    = {'Name of covariate'};
cname.strtype = 's';
cname.num     = [1 Inf];
% ---------------------------------------------------------------------
% iCC Centering
% ---------------------------------------------------------------------
iCC         = cfg_menu;
iCC.tag     = 'iCC';
iCC.name    = 'Centering';
iCC.help    = {''};
iCC.labels = {
    'Overall mean'
    'No centering'
    }';
iCC.values = {1 5};
iCC.val    = {1};
% ---------------------------------------------------------------------
% mcov Covariate
% ---------------------------------------------------------------------
mcov         = cfg_branch;
mcov.tag     = 'mcov';
mcov.name    = 'Covariate';
mcov.val     = {c cname iCC };
mcov.help    = {'Add a new covariate to your experimental design'};
% ---------------------------------------------------------------------
% generic Covariates
% ---------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Covariates';
generic.help    = {'Covariates'};
generic.values  = {mcov };
generic.num     = [0 Inf];
% ---------------------------------------------------------------------
% incint Intercept
% ---------------------------------------------------------------------
incint = cfg_menu;
incint.tag = 'incint';
incint.name = 'Intercept';
incint.help = {['By default, an intercept is always added to the model. If the ',...
    'covariates supplied by the user include a constant effect, the ',...
    'intercept may be omitted.']};
incint.labels = {'Include Intercept','Omit Intercept'};
incint.values = {1,0};
incint.val    = {1};
% ---------------------------------------------------------------------
% mreg Multiple regression
% ---------------------------------------------------------------------
mreg         = cfg_branch;
mreg.tag     = 'mreg';
mreg.name    = 'Multiple regression';
mreg.val     = {scans generic incint};
mreg.help    = {''};
% ---------------------------------------------------------------------
% name Name
% ---------------------------------------------------------------------
name         = cfg_entry;
name.tag     = 'name';
name.name    = 'Name';
name.help    = {'Name of factor, eg. ''Repetition'' '};
name.strtype = 's';
name.num     = [1 Inf];
% ---------------------------------------------------------------------
% levels Levels
% ---------------------------------------------------------------------
levels         = cfg_entry;
levels.tag     = 'levels';
levels.name    = 'Levels';
levels.help    = {'Enter number of levels for this factor, eg. 2'};
levels.strtype = 'e';
levels.num     = [Inf 1];
% ---------------------------------------------------------------------
% fact Factor
% ---------------------------------------------------------------------
fact         = cfg_branch;
fact.tag     = 'fact';
fact.name    = 'Factor';
fact.val     = {name levels dept variance gmsca ancova };
fact.help    = {'Add a new factor to your experimental design'};

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
contrast_p_value.name    = 'Contrast uncorrected p_value';
contrast_p_value.tag     = 'contrast_p_value';
contrast_p_value.strtype = 'r';
contrast_p_value.num     = [1 1];
contrast_p_value.val     = {0.05};
contrast_p_value.help    = {'Contrast uncorrected p_value'};


contrast_figures      = cfg_menu;
contrast_figures.tag  = 'contrast_figures';
contrast_figures.name = 'Generate figures';
contrast_figures.labels = {'No','Both .fig and .tiff','Only .fig','Only .tiff'};
contrast_figures.values = {0,1,2,3};
contrast_figures.val = {3};
contrast_figures.help = {'Generate contrast figures. '
    'Note .fig colorbar is incorrect - it is not saved properly by Matlab.'
    'Use .tiff to view colorbar for t-stat.'}';

figures_visible      = cfg_menu;
figures_visible.tag  = 'figures_visible';
figures_visible.name = 'Make figures visible';
figures_visible.labels = {'Yes','No'};
figures_visible.values = {1,0};
figures_visible.val = {0};
figures_visible.help = {'Make figures visible during processing.'}';

SmallFigures      = cfg_menu;
SmallFigures.tag  = 'SmallFigures';
SmallFigures.name = 'Large or small figures';
SmallFigures.labels = {'Large','Small'};
SmallFigures.values = {0,1};
SmallFigures.val = {1};
SmallFigures.help = {'Write to disk large or small (compressed) figures.'}';

colorbar_max         = cfg_entry;
colorbar_max.name    = 'Colorbar maximum value';
colorbar_max.tag     = 'colorbar_max';
colorbar_max.strtype = 'r';
colorbar_max.num     = [1 1];
colorbar_max.val     = {5};
colorbar_max.help    = {'Enter maximum value for colorbar'};

colorbar_min         = cfg_entry;
colorbar_min.name    = 'Colorbar minimum value';
colorbar_min.tag     = 'colorbar_min';
colorbar_min.strtype = 'r';
colorbar_min.num     = [1 1];
colorbar_min.val     = {2};
colorbar_min.help    = {'Enter minimum value for colorbar'};


colorbar_max2         = cfg_entry;
colorbar_max2.name    = 'Colorbar maximum value';
colorbar_max2.tag     = 'colorbar_max2';
colorbar_max2.strtype = 'r';
colorbar_max2.num     = [1 1];
colorbar_max2.val     = {-2};
colorbar_max2.help    = {'Enter maximum value for colorbar for negative maps'};

colorbar_min2         = cfg_entry;
colorbar_min2.name    = 'Colorbar minimum value';
colorbar_min2.tag     = 'colorbar_min2';
colorbar_min2.strtype = 'r';
colorbar_min2.num     = [1 1];
colorbar_min2.val     = {-5};
colorbar_min2.help    = {'Enter minimum value for colorbar for negative maps'};

colorbar_override      = cfg_branch;
colorbar_override.name      = 'Override colorbar';
colorbar_override.tag       = 'colorbar_override';
colorbar_override.val       = {colorbar_min colorbar_max colorbar_min2 colorbar_max2};
colorbar_override.help      = {'Override colorbar.'};

colorbar_default      = cfg_branch;
colorbar_default.name      = 'Default colorbar';
colorbar_default.tag       = 'colorbar_default';
colorbar_default.val       = {};
colorbar_default.help      = {'Default colorbar.'};

override_colorbar           = cfg_choice;
override_colorbar.name      = 'Override colorbar';
override_colorbar.tag       = 'override_colorbar';
override_colorbar.values    = {colorbar_default colorbar_override};
override_colorbar.val       = {colorbar_default};
override_colorbar.help      = {'Override default treatment of colorbar.'
    'User can then specify maximum and minimum values for the colorbar.'}';

GenerateInverted      = cfg_menu;
GenerateInverted.tag  = 'GenerateInverted';
GenerateInverted.name = 'Generate Inverted Responses';
GenerateInverted.labels = {'Yes','No'};
GenerateInverted.values = {1,0};
GenerateInverted.val = {1};
GenerateInverted.help = {'Generate contrasts for inverted responses.'};

GroupColorbars      = cfg_menu;
GroupColorbars.tag  = 'GroupColorbars';
GroupColorbars.name = 'Group or separate colorbars';
GroupColorbars.labels = {'Group','Separate'};
GroupColorbars.values = {1,0};
GroupColorbars.val = {0};
GroupColorbars.help = {'When considering deactivations (inverted responses),'
    'This allows displaying one (group) or two (separate) colorbars'}';

% Executable Branch
liom_1way_anova      = cfg_exbranch;
liom_1way_anova.name = 'Liom 1-way Anova Estimation';
liom_1way_anova.tag  = 'liom_1way_anova';
liom_1way_anova.val  = {NIRSmat anova_dir_name anova_level level_repeat contrast_figures contrast_p_value ...
    GroupColorbars override_colorbar figures_visible SmallFigures}; % factorial_design};
liom_1way_anova.prog = @nirs_run_liom_1way_anova;
liom_1way_anova.vout = @nirs_cfg_vout_liom_1way_anova;
liom_1way_anova.help = {'Liom 1way anova estimation.'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vout = nirs_cfg_vout_liom_1way_anova(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});