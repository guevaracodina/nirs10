function wls_bglm_specify = nirs_run_liom_GLM_specify_cfg
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'Stat');
target_sampling_rate = nirs_dfg_target_sampling_rate;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vasomotion
%%%%%%%%%%%%%%%%%%%%%%%%%%%
select_chromophore         = cfg_menu;
select_chromophore.tag     = 'select_chromophore';
select_chromophore.name    = 'Select chromophore';
select_chromophore.help    = {'Select chromophore'};
select_chromophore.labels = {'HbT' 'HbR' 'HbO' 'HbO&HbR'};
select_chromophore.values  = {  1   2   3   4};
select_chromophore.val = {3};

vasomotion_on         = cfg_branch;
vasomotion_on.tag     = 'vasomotion_on';
vasomotion_on.name    = 'Include vasomotion regressor';
vasomotion_on.val     = {select_chromophore};
vasomotion_on.help    = {'Include vasomotion regressor.'};

no_vasomotion         = cfg_branch;
no_vasomotion.tag     = 'no_vasomotion';
no_vasomotion.name    = 'No vasomotion regressor';
no_vasomotion.val     = {};
no_vasomotion.help    = {};

vasomotion_choice        = cfg_choice;
vasomotion_choice.name   = 'Choose vasomotion method';
vasomotion_choice.tag    = 'vasomotion_choice';
vasomotion_choice.values = {no_vasomotion,vasomotion_on};
vasomotion_choice.val    = {no_vasomotion};
vasomotion_choice.help   = {'Choose whether to include a vasomotion regressor in the GLM.'}';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General Linear Model Specification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------------
% dir Directory
% ---------------------------------------------------------------------
dir         = cfg_files;
dir.tag     = 'dir';
dir.name    = 'Stats Directory';
dir.help    = {'Select a directory where the NIRS_SPM.mat files containing the specified design matrix will be written.'};
dir.filter = 'dir';
dir.ufilter = '.*';
dir.num     = [1 1];

% ---------------------------------------------------------------------
% sessions Sessions
% ---------------------------------------------------------------------
sessions         = cfg_entry;
sessions.name    = 'Sessions';
sessions.tag     = 'sessions';
sessions.strtype = 'r';
sessions.num     = [0 Inf];
sessions.val     = {''};
sessions.help    = {'Enter the numbers of the sessions to include. If not specified, all sessions included in the NIRS matrix will be included.'};

% ---------------------------------------------------------------------
% units Units for design
% ---------------------------------------------------------------------
units         = cfg_menu;
units.tag     = 'units';
units.name    = 'Units for design';
units.help    = {'The onsets of events or blocks can be specified in either scans or seconds.'};
units.labels = {
    'Scans'
    'Seconds'
    };
units.values  = {
    0
    1
    };
units.val = {1};

time_res      = cfg_entry;
time_res.tag  = 'time_res';
time_res.name = 'Time resolution';
time_res.val = {1};
time_res.strtype = 'r';
time_res.num     = [1 1];
time_res.help    = {'Time resolution for onsets will be given by NIRS sampling rate divided by this factor  - value is 10 in NIRS_SPM.'};

input_onsets         = cfg_files;
input_onsets.name    = 'Select onset files'; % The displayed name
input_onsets.tag     = 'input_onsets';
input_onsets.filter  = 'mat';
input_onsets.val{1}  = {''};
input_onsets.num     = [0 Inf];     % Number of inputs required
input_onsets.help    = {'Select onset files for each session of this subject.'}; % help text displayed

% ---------------------------------------------------------------------
% multi_reg Multiple regressors
% ---------------------------------------------------------------------
multi_reg         = cfg_files;
multi_reg.tag     = 'multi_reg';
multi_reg.name    = 'Multiple regressors';
multi_reg.val{1} = {''};
multi_reg.help    = {
    'Select the *.mat/*.txt file containing details of your multiple regressors. '
    ''
    'If you have multiple regressors eg. realignment parameters, then entering the details a regressor at a time is very inefficient. This option can be used to load all the required information in one go. '
    ''
    'You will first need to create a *.mat file containing a matrix R or a *.txt file containing the regressors. Each column of R will contain a different regressor. When SPM creates the design matrix the regressors will be named R1, R2, R3, ..etc.'
    }';
multi_reg.filter = 'mat';
multi_reg.ufilter = '.*';
multi_reg.num     = [0 Inf];


subj         = cfg_branch;
subj.tag     = 'subj';
subj.name    = 'Subject';
subj.val     = {input_onsets multi_reg};
subj.help    = {};



% ---------------------------------------------------------------------
% Noise method
% ---------------------------------------------------------------------
nirs_noise         = cfg_menu;
nirs_noise.tag     = 'nirs_noise';
nirs_noise.name    = 'Noise method';
nirs_noise.help    = {'Choose method for noise treatment.'}';
nirs_noise.labels  = {
    'precoloring'
    'prewhitening'
    };
nirs_noise.values  = {
    0
    1
    };
nirs_noise.def = @(val)nirs_get_defaults(...
    'model_specify.wls_bglm_specify.wls_or_bglm.NIRS_SPM.nirs_noise', val{:});

%Get bases for HRF
bases = nirs_dfg_hrf;

% ---------------------------------------------------------------------
% volt Model Interactions (Volterra)
% ---------------------------------------------------------------------
volt         = cfg_menu;
volt.tag     = 'volt';
volt.name    = 'Model Interactions (Volterra)';
volt.help    = {
    'Generalized convolution of inputs (U) with basis set (bf).'
    ''
    'For first order expansions the causes are simply convolved (e.g. stick functions) in U.u by the basis functions in bf to create a design matrix X.  For second order expansions new entries appear in ind, bf and name that correspond to the interaction among the orginal causes. The basis functions for these efects are two dimensional and are used to assemble the second order kernel. Second order effects are computed for only the first column of U.u.'
    'Interactions or response modulations can enter at two levels.  Firstly the stick function itself can be modulated by some parametric variate (this can be time or some trial-specific variate like reaction time) modeling the interaction between the trial and the variate or, secondly interactions among the trials themselves can be modeled using a Volterra series formulation that accommodates interactions over time (and therefore within and between trial types).'
    }';
volt.labels = {
    'Do not model Interactions'
    'Model Interactions (2nd Volterra)'
    'Model 3rd Volterra'
    }';
volt.values = {1 2 3};
volt.val = {1};

% ----------------------------------------------------------------------
% Early response
% ----------------------------------------------------------------------

time_shift         = cfg_entry;
time_shift.tag     = 'time_shift';
time_shift.name    = 'Time shift of onsets (early response)';
time_shift.strtype = 'r';
time_shift.help    = {
    'This option is for the study of the early response of onsets.'
    'Please specify the desired shifting time (in seconds) here (e.g. 0, 3, 6, 9, 12, 15). '
    'If a nonzero real value is entered, the timing of the onsets of interest will be adjusted ahead.'
    }';
time_shift.num     = [1 1];
time_shift.val     = {0};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LIOM General Linear Model Specification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GenerateHbT      = cfg_menu;
GenerateHbT.tag  = 'GenerateHbT';
GenerateHbT.name = 'Generate HbT';
GenerateHbT.labels = {'Yes','No'};
GenerateHbT.values = {1,0};
GenerateHbT.def = @(val)nirs_get_defaults(...
    'model_specify.wls_bglm_specify.GenerateHbT', val{:});
GenerateHbT.help = {'Generate HbT.'};

flag_window      = cfg_menu;
flag_window.tag  = 'flag_window';
flag_window.name = 'Show Design Matrix';
flag_window.labels = {'Yes','No'};
flag_window.values = {1,0};
flag_window.def = @(val)nirs_get_defaults(...
    'model_specify.wls_bglm_specify.flag_window', val{:});
flag_window.help = {'Show design matrix.'};

WLS_J0         = cfg_entry;
WLS_J0.name    = 'Wavelet depth J0';
WLS_J0.tag     = 'WLS_J0';
WLS_J0.strtype = 'r';
WLS_J0.num     = [1 1];
WLS_J0.def     = @(val)nirs_get_defaults(...
    'model_specify.wls_bglm_specify.wls_or_bglm.WLS.WLS_J0', val{:});
WLS_J0.help    = {'Enter wavelet depth J0.'};

WLS_L0         = cfg_entry;
WLS_L0.name    = 'Wavelet depth L0';
WLS_L0.tag     = 'WLS_L0';
WLS_L0.strtype = 'r';
WLS_L0.num     = [1 1];
WLS_L0.def     = @(val)nirs_get_defaults(...
    'model_specify.wls_bglm_specify.wls_or_bglm.WLS.WLS_L0', val{:});
WLS_L0.help    = {'Enter wavelet depth L0.'};

WLS_threshold_drift         = cfg_entry;
WLS_threshold_drift.name    = 'Wavelet correlation threshold for drifts';
WLS_threshold_drift.tag     = 'WLS_threshold_drift';
WLS_threshold_drift.strtype = 'r';
WLS_threshold_drift.num     = [1 1];
WLS_threshold_drift.def     = @(val)nirs_get_defaults(...
    'model_specify.wls_bglm_specify.wls_or_bglm.WLS.WLS_threshold_drift', val{:});
WLS_threshold_drift.help    = {'Enter wavelet correlation threshold for drifts.'};

WLS         = cfg_branch;
WLS.tag     = 'WLS';
WLS.name    = 'Wavelet least-squares';
WLS.val     = {WLS_J0 WLS_threshold_drift WLS_L0};
WLS.help    = {'Specify options for wavelet least-squares method.'};

BGLM_fmax         = cfg_entry;
BGLM_fmax.name    = 'Maximum frequency for drifts';
BGLM_fmax.tag     = 'BGLM_fmax';
BGLM_fmax.strtype = 'r';
BGLM_fmax.num     = [1 1];
BGLM_fmax.def     = @(val)nirs_get_defaults(...
    'model_specify.wls_bglm_specify.wls_or_bglm.BGLM.BGLM_fmax', val{:});
BGLM_fmax.help    = {'Enter maximum frequency for drifts in Hz.'};

BGLM_degre         = cfg_entry;
BGLM_degre.name    = 'Polynomial degree for drifts';
BGLM_degre.tag     = 'BGLM_degre';
BGLM_degre.strtype = 'r';
BGLM_degre.num     = [1 1];
BGLM_degre.def     = @(val)nirs_get_defaults(...
    'model_specify.wls_bglm_specify.wls_or_bglm.BGLM.BGLM_degre', val{:});
BGLM_degre.help    = {'Enter polynomial degree for drifts.'};

BGLM_threshold_drift         = cfg_entry;
BGLM_threshold_drift.name    = 'Threshold for drifts';
BGLM_threshold_drift.tag     = 'BGLM_threshold_drift';
BGLM_threshold_drift.strtype = 'r';
BGLM_threshold_drift.num     = [1 1];
BGLM_threshold_drift.def     = @(val)nirs_get_defaults(...
    'model_specify.wls_bglm_specify.wls_or_bglm.BGLM.BGLM_threshold_drift', val{:});
BGLM_threshold_drift.help    = {'Enter correlation threshold for drifts.'};

BGLM         = cfg_branch;
BGLM.tag     = 'BGLM';
BGLM.name    = 'Bayesian GLM';
BGLM.val     = {BGLM_fmax BGLM_degre BGLM_threshold_drift};
BGLM.help    = {'Specify options for Bayesian GLM method.'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.6 High pass and low pass filters (optional)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------------------------------------------------------------------
% hpf High-pass filter
% ---------------------------------------------------------------------

hpf_wavelet_iter      = cfg_entry;
hpf_wavelet_iter.tag  = 'hpf_wavelet_iter';
hpf_wavelet_iter.name = 'Wavelet iterations';
hpf_wavelet_iter.val = {4};
%hpf_wavelet_iter.def    = @(val)nirs_get_defaults('hpf_wavelet_iter', val{:});
hpf_wavelet_iter.strtype = 'r';
hpf_wavelet_iter.num     = [1 1];
hpf_wavelet_iter.help    = {'Specify wavelet iterations - default is 4 in NIRS_SPM.'};

hpf_wavelet = cfg_branch;
hpf_wavelet.tag     = 'hpf_wavelet';
hpf_wavelet.name    = 'Wavelet Filter';
hpf_wavelet.val     = {hpf_wavelet_iter};
hpf_wavelet.help    = {'Specify properties of wavelet filter'};

hpf_dct_cutoff      = cfg_entry;
hpf_dct_cutoff.tag  = 'hpf_dct_cutoff';
hpf_dct_cutoff.name = 'DCT cutoff in seconds';
hpf_dct_cutoff.val  = {128};
hpf_dct_cutoff.strtype = 'r';
hpf_dct_cutoff.num     = [1 1];
hpf_dct_cutoff.help    = {'Specify DCT cutoff in seconds.'};

hpf_dct = cfg_branch;
hpf_dct.tag     = 'hpf_dct';
hpf_dct.name    = 'DCT Filter';
hpf_dct.val     = {hpf_dct_cutoff};
hpf_dct.help    = {'Specify properties of Discrete Cosine Transform filter'};

hpf_none         = cfg_branch;
hpf_none.tag     = 'hpf_none';
hpf_none.name    = 'No high pass filter';
hpf_none.help    = {'No high pass filter.'};

nirs_hpf           = cfg_choice;
nirs_hpf.name      = 'High-pass filter';
nirs_hpf.tag       = 'nirs_hpf';
nirs_hpf.values    = {hpf_none
    hpf_wavelet
    hpf_dct};
nirs_hpf.val       = {hpf_none};
nirs_hpf.help      = {'Choose high-pass filter.'};

% ---------------------------------------------------------------------
% lpf Low-pass filter
% ---------------------------------------------------------------------
fwhm1      = cfg_entry;
fwhm1.tag  = 'fwhm1';
fwhm1.name = 'FWHM in seconds';
fwhm1.val = {1.5};
%fwhm1.def    = @(val)nirs_get_defaults('fwhm1', val{:});
fwhm1.strtype = 'r';
fwhm1.num     = [1 1];
%fwhm1.def = @(val)nirs_get_defaults('configMC1.scalpPpties_l2', val{:});
fwhm1.help    = {'FWHM in seconds.'};

lpf_gauss         = cfg_branch;
lpf_gauss.tag     = 'lpf_gauss';
lpf_gauss.name    = 'Gaussian Filter';
lpf_gauss.val     = {fwhm1};
lpf_gauss.help    = {'Specify properties of Gaussian filter'};

lpf_none         = cfg_branch;
lpf_none.tag     = 'lpf_none';
lpf_none.name    = 'No low pass filter';
lpf_none.help    = {'No low pass filter.'};

lpf_hrf         = cfg_branch;
lpf_hrf.tag     = 'lpf_hrf';
lpf_hrf.name    = 'HRF Filter';
lpf_hrf.help    = {'HRF filter'};

nirs_lpf           = cfg_choice;
nirs_lpf.name      = 'Low-pass filter';
nirs_lpf.tag       = 'nirs_lpf';
nirs_lpf.values    = {lpf_none
    lpf_gauss
    lpf_hrf};
nirs_lpf.val       = {lpf_hrf};
nirs_lpf.help      = {'Choose low-pass filter.'
    'Currently, the only valid choice is the LPF based on the hemodynamic response function.'
    '(Or a Gaussian filter with a long cutoff of 4 or 5 s.)'
    'A Gaussian filter with a short cutoff leads to an inaccurate estimation '
    'of the number of degrees of freedom and hence incorrect statistics (t-stats being larger than they should be)'}';

NIRS_SPM         = cfg_branch;
NIRS_SPM.tag     = 'NIRS_SPM';
NIRS_SPM.name    = 'NIRS_SPM MDL';
NIRS_SPM.val     = {nirs_noise nirs_hpf nirs_lpf};
NIRS_SPM.help    = {'Specify options for NIRS_SPM minimum description length(MDL).'};

wls_or_bglm      = cfg_choice;
wls_or_bglm.tag  = 'wls_or_bglm';
wls_or_bglm.name = 'WLS, BGLM, NIRS_SPM';
%wls_or_bglm.labels = {'WLS','BGLM', 'NIRS_SPM'};
%wls_or_bglm.values = {1,2,3};
wls_or_bglm.values = {WLS,BGLM,NIRS_SPM};
wls_or_bglm.val  = {NIRS_SPM};
%wls_or_bglm.def = @(val)nirs_get_defaults('model_specify.wls_bglm_specify.wls_or_bglm', val{:});
wls_or_bglm.help = {'Choose which GLM method to use:'
    'WLS: wavelet least square'
    'BGLM: Bayesian general linear model'
    'NIRS_SPM: Ye et al methods (MDL), with either precoloring or prewhitening.'}';

GLM_include_cardiac    = cfg_menu;
GLM_include_cardiac.name   = 'Include cardiac regressor';
GLM_include_cardiac.tag    = 'GLM_include_cardiac';
GLM_include_cardiac.labels = {'Yes','No'};
GLM_include_cardiac.values = {1,0};
GLM_include_cardiac.def    = @(val)nirs_get_defaults(...
    'model_specify.wls_bglm_specify.GLM_include_cardiac', val{:});
GLM_include_cardiac.help   = {'Include cardiac regressor if available.'};

GLM_include_Mayer    = cfg_menu;
GLM_include_Mayer.name   = 'Include Mayer wave regressor';
GLM_include_Mayer.tag    = 'GLM_include_Mayer';
GLM_include_Mayer.labels = {'Yes','No'};
GLM_include_Mayer.values = {1,0};
GLM_include_Mayer.def    = @(val)nirs_get_defaults(...
    'model_specify.wls_bglm_specify.GLM_include_Mayer', val{:});
GLM_include_Mayer.help   = {'Include Mayer wave regressor if available.'};

channel_pca      = cfg_menu;
channel_pca.tag  = 'channel_pca';
channel_pca.name = 'Spatial Principal Component Removal';
channel_pca.labels = {'Yes','No'};
channel_pca.values = {1,0};
channel_pca.def = @(val)nirs_get_defaults('model_specify.wls_bglm_specify.channel_pca', val{:});
channel_pca.help = {'Choose whether to do a channel PCA removal: '
    'Principal component analysis and removing the largest eigenvalue.'}';

NumPCAComponents         = cfg_entry;
NumPCAComponents.name    = 'Number of PCA components to remove';
NumPCAComponents.tag     = 'NumPCAComponents';
NumPCAComponents.strtype = 'r';
NumPCAComponents.num     = [1 1];
NumPCAComponents.val     = {1};
NumPCAComponents.help    = {'Enter number of PCA components to be removed.'
    'This option will only be used when the PCA option above is selected'}';

hpf_butter_freq         = cfg_entry;
hpf_butter_freq.name    = 'Cutoff frequency for HPF';
hpf_butter_freq.tag     = 'hpf_butter_freq';
hpf_butter_freq.strtype = 'r';
hpf_butter_freq.num     = [1 1];
hpf_butter_freq.val     = {0.01};
hpf_butter_freq.help    = {'Enter cutoff frequency in Hz for Butterworth HPF.'};

hpf_butter_order         = cfg_entry;
hpf_butter_order.name    = 'Order of Butterworth HPF';
hpf_butter_order.tag     = 'hpf_butter_order';
hpf_butter_order.strtype = 'r';
hpf_butter_order.num     = [1 1];
hpf_butter_order.val     = {2};
hpf_butter_order.help    = {'Enter order of Butterworth HPF (preferred value = 2).'};

hpf_butter_On         = cfg_branch;
hpf_butter_On.tag     = 'hpf_butter_On';
hpf_butter_On.name    = 'Butterworth HP filter';
hpf_butter_On.val     = {hpf_butter_freq hpf_butter_order};
hpf_butter_On.help    = {'Butterworth high-pass filter.'};

hpf_butter_Off         = cfg_branch;
hpf_butter_Off.tag     = 'hpf_butter_Off';
hpf_butter_Off.name    = 'HP filter off';
hpf_butter_Off.val     = {};
hpf_butter_Off.help    = {'High pass filter turned off.'};

remove_linear         = cfg_branch;
remove_linear.tag     = 'remove_linear';
remove_linear.name    = 'Remove linear trend only';
remove_linear.val     = {};
remove_linear.help    = {'Remove linear trend only.'};

GLM_remove_linear         = cfg_branch;
GLM_remove_linear.tag     = 'GLM_remove_linear';
GLM_remove_linear.name    = 'Remove linear trend using the GLM';
GLM_remove_linear.val     = {};
GLM_remove_linear.help    = {'Remove linear trend using the GLM.'
    'That is, a regressor will be added to the GLM for the linear trend.'}';

SPM_cosine_filter         = cfg_branch;
SPM_cosine_filter.tag     = 'SPM_cosine_filter';
SPM_cosine_filter.name    = 'Remove trend using cosines in the GLM';
SPM_cosine_filter.val     = {};
SPM_cosine_filter.help    = {'Remove trends using cosines as done in SPM'}';

NoNIRSconfounds         = cfg_branch;
NoNIRSconfounds.tag     = 'NoNIRSconfounds';
NoNIRSconfounds.name    = 'No NIRS channels as confounds';
NoNIRSconfounds.val     = {};
NoNIRSconfounds.help    = {'No NIRS channels as confounds.'};

% NumChConfounds         = cfg_entry;
% NumChConfounds.name    = 'Maximum Number of Confounds';
% NumChConfounds.tag     = 'NumChConfounds';
% NumChConfounds.strtype = 'r';
% NumChConfounds.num     = [1 1];
% NumChConfounds.val     = {1};
% NumChConfounds.help    = {'Enter maximum number of NIRS channels to be included as physiological confounds.'};

MinChDist         = cfg_entry;
MinChDist.name    = 'Minimum Channel Length';
MinChDist.tag     = 'MinChDist';
MinChDist.strtype = 'r';
MinChDist.num     = [1 1];
MinChDist.val     = {1.5};
MinChDist.help    = {'Enter minimum channel length allowed for inclusion as confound.'};

MaxChDist         = cfg_entry;
MaxChDist.name    = 'maximum Channel Length';
MaxChDist.tag     = 'MaxChDist';
MaxChDist.strtype = 'r';
MaxChDist.num     = [1 1];
MaxChDist.val     = {2.5};
MaxChDist.help    = {'Enter maximum channel length allowed for inclusion as confound.'};

NumChConfounds         = cfg_entry;
NumChConfounds.name    = 'Maximum Number of Confounds';
NumChConfounds.tag     = 'NumChConfounds';
NumChConfounds.strtype = 'r';
NumChConfounds.num     = [1 1];
NumChConfounds.val     = {1};
NumChConfounds.help    = {'Enter maximum number of NIRS channels to be included as physiological confounds.'};

NIRSconfounds         = cfg_branch;
NIRSconfounds.tag     = 'NIRSconfounds';
NIRSconfounds.name    = 'NIRS channels as confounds';
NIRSconfounds.val     = {NumChConfounds MinChDist MaxChDist};
NIRSconfounds.help    = {'NIRS channels as confounds.'};

NIRSchannelsConfound         = cfg_choice;
NIRSchannelsConfound.tag     = 'NIRSchannelsConfound';
NIRSchannelsConfound.name    = 'NIRS channels as confounds';
NIRSchannelsConfound.values  = {NoNIRSconfounds NIRSconfounds};
NIRSchannelsConfound.val     = {NoNIRSconfounds};
NIRSchannelsConfound.help    = {'NIRS channels  as confound regressors.'
    'When using this option, selected channels will be filtered with the'
    'Same parameters as for the GLM and included as confound regressors'
    'This may be useful to remove physiological noise'
    'Only HbO channels will be used'}';

generate_trRV      = cfg_menu;
generate_trRV.tag  = 'generate_trRV';
generate_trRV.name = 'Generate TrRV';
generate_trRV.labels = {'Yes','No'};
generate_trRV.values = {1,0};
generate_trRV.val = {1};
generate_trRV.help = {'Generate TrRV and TrRVRV - needed for interpolated maps.'
    'Careful! Note that TrRV is required for the NIRS_SPM method, to calculate t-stats.'}';

TrRVRVexact      = cfg_menu;
TrRVRVexact.tag  = 'TrRVRVexact';
TrRVRVexact.name = 'Exact calculation of TrRVRV';
TrRVRVexact.labels = {'Yes','No'};
TrRVRVexact.values = {1,0};
TrRVRVexact.val = {0};
TrRVRVexact.help = {'Perform an exact calculation for TrRVRV (long for long data files)'
    'Note that discrepancies of 30% have been found between this approximation and the exact result.'
    'This will affect the number of degrees of freedom, and therefore the statistical thresholds'}';

hpf_butter      = cfg_choice;
hpf_butter.tag  = 'hpf_butter';
hpf_butter.name = 'Additional High Pass Filter';
hpf_butter.values = {hpf_butter_On hpf_butter_Off remove_linear GLM_remove_linear SPM_cosine_filter};
hpf_butter.val = {hpf_butter_On};
hpf_butter.help = {'Additional High Pass Filter'}';

%This is a duplication -- already specified in HRF -- but its purpose is
%to be able to add derivatives of a gamma function.
derivs         = cfg_menu;
derivs.tag     = 'derivs';
derivs.name    = 'Model derivatives';
derivs.help    = {'Model HRF Derivatives. The canonical HRF combined with time and dispersion derivatives comprise an ''informed'' basis set, as the shape of the canonical response conforms to the hemodynamic response that is commonly observed. The incorporation of the derivate terms allow for variations in subject-to-subject and voxel-to-voxel responses. The time derivative allows the peak response to vary by plus or minus a second and the dispersion derivative allows the width of the response to vary. The informed basis set requires an SPM{F} for inference. T-contrasts over just the canonical are perfectly valid but assume constant delay/dispersion. The informed basis set compares favourably with eg. FIR bases on many data sets. '};
derivs.labels = {
    'No derivatives'
    'Time derivatives'
    'Time and Dispersion derivatives'
    }';
derivs.values = {[0 0] [1 0] [1 1]};
derivs.val = {[0 0]};

% Executable Branch
wls_bglm_specify      = cfg_exbranch;
wls_bglm_specify.name = 'LIOM GLM Specification';
wls_bglm_specify.tag  = 'wls_bglm_specify';
wls_bglm_specify.val  = {NIRSmat redo1 NIRSmatCopyChoice sessions subj units time_res derivs bases ...
    volt time_shift GLM_include_cardiac GLM_include_Mayer vasomotion_choice NIRSchannelsConfound GenerateHbT flag_window ...
    channel_pca NumPCAComponents hpf_butter generate_trRV TrRVRVexact ... %filter_design_matrix ...
    target_sampling_rate wls_or_bglm };
wls_bglm_specify.prog = @nirs_run_liom_GLM_specify;
wls_bglm_specify.vout = @nirs_cfg_vout_liom_GLM_specify;
wls_bglm_specify.help = {'Specify LIOM General Linear Model.'};

function vout = nirs_cfg_vout_liom_GLM_specify(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});