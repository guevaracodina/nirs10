function ReMLreconstruct1 = nirs_run_ReMLreconstruct_cfg
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'ReML');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reconstructions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


temp_pts         = cfg_entry;
temp_pts.name    = 'Point temporel de l''inversion'; % The displayed name
temp_pts.tag     = 'temp_pts';       %file names
temp_pts.strtype = 'r';
temp_pts.num     = [1 Inf];     % Number of inputs required
temp_pts.val     = {1};
% temp_pts.def = @(val)nirs_get_defaults('configMC1.nphotons', val{:});
temp_pts.help    = {'Input time points.'};

specific_points           = cfg_branch;
specific_points.name      = 'Select specific_points';
specific_points.tag       = 'specific_points';
specific_points.val       = {temp_pts};
specific_points.help      = {'Select specific points.'};

downsample_freq         = cfg_entry;
downsample_freq.name    = 'Downsampling frequency in Hz'; % The displayed name
downsample_freq.tag     = 'downsample_freq';       %file names
downsample_freq.strtype = 'r';
downsample_freq.num     = [1 1];     % Number of inputs required
downsample_freq.val     = {1};
downsample_freq.help    = {'Enter downsampling frequency in Hz.'
    'A target sampling frequency will be generated, which may however'
    'be only approximately equal to the specified downsampling frequency,'
    'but it will correspond to the actual frequency of selecting every Nth point'}';

all_points_downsampled           = cfg_branch;
all_points_downsampled.name      = 'Select all points, downsampled';
all_points_downsampled.tag       = 'all_points_downsampled';
all_points_downsampled.val       = {downsample_freq};
all_points_downsampled.help      = {'Select a subset of all points,'
    'downsampled to a set frequency'}';

psel_choice           = cfg_choice;
psel_choice.name      = 'Temporal point selection method';
psel_choice.tag       = 'psel_choice';
psel_choice.values    = {all_points_downsampled specific_points};
psel_choice.val       = {all_points_downsampled};
psel_choice.help      = {'Temporal point selection method'}';

sens_vxsize= cfg_entry;
sens_vxsize.name    = 'Voxel size in sensitivity matrix'; % The displayed name
sens_vxsize.tag     = 'sens_vxsize';       %file names
sens_vxsize.strtype = 'r';
sens_vxsize.num     = [1 1];     % Number of inputs required
sens_vxsize.val     = {5};
% sens_vxsize.def = @(val)nirs_get_defaults('configMC1.nphotons', val{:});
sens_vxsize.help    = {'Enter voxel size for reconstruction, in millimeters.'};

anat_segT1         = cfg_files; %Select MC segmented volume for this subject
anat_segT1.name    = 'Anatomical segmented image'; % The displayed name
anat_segT1.tag     = 'anat_segT1';       %file names
anat_segT1.filter  = 'image';
anat_segT1.ufilter = '.nii';
anat_segT1.num     = [1 1];     % Number of inputs required
anat_segT1.help    = {'Anatomical segmented image with NewSegment and MCsegment.'}; % help text displayed

% Priors
hb_relative_evolution        = cfg_menu;
hb_relative_evolution.name   = 'Relative evolution HbO/HbR';
hb_relative_evolution.tag    = 'hb_relative_evolution';
hb_relative_evolution.labels = {'Yes','No'};
hb_relative_evolution.values = {1,0};
%hb_relative_evolution.def    = @(val)nirs_get_defaults('hb_relative_evolution', val{:});
hb_relative_evolution.help   = {'Choose type of configuration files to generate.'};

priors      = cfg_branch;
priors.name = 'Priors';
priors.tag  = 'priors';
priors.val  = {anat_segT1 hb_relative_evolution};
priors.help = {'Choose priors you want to use for the reconstruction.'};

WLruns           = cfg_menu;
WLruns.name      = 'Runs';
WLruns.tag       = 'WLruns';
WLruns.labels    = {'One' 'Each WL separately'};
WLruns.values    = {1,2};
WLruns.val       = {1};
WLruns.help      = {''};

beta_wtd           = cfg_menu;
beta_wtd.name      = 'Beta : mua or Hbs';
beta_wtd.tag       = 'beta_wtd';
beta_wtd.labels    = {'mua' 'hbs'};
beta_wtd.values    = {1,2};
beta_wtd.val       = {1};
beta_wtd.help      = {'Choose Delta mua or Delta[HbO] and Delta[HBR]'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration: 3D reconstruction -- Tikhonov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha         = cfg_entry;
alpha.name    = 'Hyperparameter';
alpha.tag     = 'alpha';
alpha.strtype = 'r';
alpha.num     = [1 1];
alpha.val     = {1};
alpha.help    = {'Tunes the model :'
    'If small, the solution favors minimizing the residual error with the measured data.'
    'If large, the solution is biased towards matching the prior (no activity).'}';

alpha2         = cfg_entry;
alpha2.name    = 'Hyperparameter Mask';
alpha2.tag     = 'alpha2';
alpha2.strtype = 'r';
alpha2.num     = [1 1];
alpha2.val     = {100};
alpha2.help    = {'Change weight for mask. Actual weight used in the code will be alpha*alpha2.'}';

wgmc      = cfg_branch;
wgmc.name = 'Constraint on White and Grey Matter';
wgmc.tag  = 'wgmc';
wgmc.help = {''};

samcs      = cfg_branch;
samcs.name = 'Same as Monte Carlo simulation';
samcs.tag  = 'samcs';
samcs.help = {'Helmet informations will be extracted from ''.nirs'' file.'};

timask         = cfg_files;
timask.name    = 'Timask';
timask.tag     = 'timask';
timask.filter  = 'image';
timask.ufilter = '.nii';
timask.num     = [1 1];
timask.help    = {'.'};

tikh_mask         = cfg_choice;
tikh_mask.tag     = 'tikh_mask';
tikh_mask.name    = 'Mask for constraint';
tikh_mask.values  = {wgmc samcs timask};
tikh_mask.val     = {wgmc};
tikh_mask.help    = {'Choose mask on witch you want to constraint reconstruction.'
    'It can be the same as the one selected for Monte Carlo simulations then choose ''same as Monte Carlo simulations'' or any image.'}';

tikh_SC      = cfg_branch;
tikh_SC.name = 'Constraints';
tikh_SC.tag  = 'tikh_SC';
tikh_SC.val  = {tikh_mask alpha2};
tikh_SC.help = {'Choose mask and hyperparameter.'};

spatial_constraint         = cfg_repeat;
spatial_constraint.tag     = 'generic1';
spatial_constraint.name    = 'Spatial constraint for extended Tikhonov inversion';
spatial_constraint.help    = {'Spatial constraint'}';
spatial_constraint.values  = {tikh_SC};
spatial_constraint.num     = [0 2];

tikhonov      = cfg_branch;
tikhonov.name = 'Tikhonov';
tikhonov.tag  = 'tikhonov';
tikhonov.help = {''};

simple_bayes      = cfg_branch;
simple_bayes.name = 'Simple Bayesian Interpretation';
simple_bayes.tag  = 'simple_bayes';
simple_bayes.help = {''};

tikh_method         = cfg_choice;
tikh_method.tag     = 'tikh_method';
tikh_method.name    = 'Tikhonov regularization method';
tikh_method.values  = {tikhonov simple_bayes};
tikh_method.val     = {tikhonov};
tikh_method.help    = {'Choose Tikhonov regularization reconstruction method (all taken from Hierarchical Bayesian regularization of reconstructions for diffuse optical tomography using multiple priors, Huppert).'
    '-- Tikhonov is the simplest method of regularization'
    '-- Extended Tikhonov uses Li et al. model'
    '-- Simple Bayesian Interpretation uses covariances for the norms.'
    }';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration: 3D reconstruction -- ReML reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ReML_method           = cfg_menu;
ReML_method.name      = 'ReML method';
ReML_method.tag       = 'ReML_method';
ReML_method.labels    = {'Huppert' 'SPM'};
ReML_method.values    = {0,1};
ReML_method.val       = {0};
ReML_method.help      = {'Choose ReML reconstruction method.'};

%%%%%%%%%%%%%BOLD
dir_in         = cfg_files;
dir_in.tag     = 'dir_in';
dir_in.name    = 'MonteCarlo output directory';
dir_in.help    = {'Select the MonteCarlo simulation output directory.'};
dir_in.filter = 'dir';
dir_in.val{1} = {''};
dir_in.ufilter = '.*';
dir_in.num     = [0 1];

% Executable Branch
ReMLreconstruct1      = cfg_exbranch;
ReMLreconstruct1.name = '3D NIRS data ReML reconstruction';
ReMLreconstruct1.tag  = 'ReMLreconstruct1';
ReMLreconstruct1.val  = {NIRSmat redo1 NIRSmatCopyChoice beta_wtd psel_choice dir_in sens_vxsize ReML_method WLruns};
ReMLreconstruct1.prog = @nirs_run_ReMLreconstruct;
ReMLreconstruct1.vout = @nirs_cfg_vout_ReMLreconstruct;
ReMLreconstruct1.help = {'Run 3D NIRS data reconstruction.'};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_ReMLreconstruct(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
