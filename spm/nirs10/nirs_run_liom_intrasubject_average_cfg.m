function liom_intrasubject_average = nirs_run_liom_intrasubject_average_cfg
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'Avg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Averaging method
%%%%%%%%%%%%%%%%%%%%%%%%%%%

onset_delays         = cfg_entry;
onset_delays.name    = 'Onset delays in seconds';
onset_delays.tag     = 'onset_delays';
onset_delays.strtype = 'r';
onset_delays.num     = [1 Inf];
onset_delays.val     = {3};
onset_delays.help    = {'Enter an array of onset delays in seconds; negative numbers are allowed.'};

onset_duration         = cfg_entry;
onset_duration.name    = 'Onset duration in seconds';
onset_duration.tag     = 'onset_duration';
onset_duration.strtype = 'r';
onset_duration.num     = [1 1];
onset_duration.val     = {2};
onset_duration.help    = {'Enter onset duration.'
    'For an event related protocol, this could be 1 or 2 seconds. '
    'For a block protocol, this will typically be the length of the block.'}';

seizure_evolution         = cfg_branch;
seizure_evolution.tag     = 'seizure_evolution';
seizure_evolution.name    = 'Seizure Evolution';
seizure_evolution.val     = {onset_delays onset_duration};
seizure_evolution.help    = {'Seizure evolution.'
    'Goal here is to use Z-scores where the standard deviation is calculated from the baseline'.'}';

onset_delay         = cfg_entry;
onset_delay.name    = 'Onset delay in seconds';
onset_delay.tag     = 'onset_delay';
onset_delay.strtype = 'r';
onset_delay.num     = [1 1];
onset_delay.val     = {3};
onset_delay.help    = {'Enter onset delay in seconds -- typically 3 seconds.'};

onset_duration         = cfg_entry;
onset_duration.name    = 'Onset duration in seconds';
onset_duration.tag     = 'onset_duration';
onset_duration.strtype = 'r';
onset_duration.num     = [1 1];
onset_duration.val     = {2};
onset_duration.help    = {'Enter onset duration.'
    'For an event related protocol, this could be 1 or 2 seconds. '
    'For a block protocol, this will typically be the length of the block.'}';

block_averaging         = cfg_branch;
block_averaging.tag     = 'block_averaging';
block_averaging.name    = 'Block averaging';
block_averaging.val     = {onset_delay onset_duration};
block_averaging.help    = {'Block averaging.'
    'This works by first averaging the data for each onset over the specified window'
    'These averages for each onset are then averaged and their standard deviation is calculated.'}';

SL_List         = cfg_entry;
SL_List.name    = 'List of sessions to average';
SL_List.tag     = 'SL_List';
SL_List.strtype = 'r';
SL_List.num     = [0 Inf];
SL_List.val     = {''};
SL_List.help    = {'Enter session numbers of sessions to be averaged.'}';

SL_Label         = cfg_entry; 
SL_Label.name    = 'Name for this session list';
SL_Label.tag     = 'SL_Label';
SL_Label.strtype = 's';
SL_Label.num     = [0 Inf];
SL_Label.val{1}  = '';
SL_Label.help    = {'Enter name or label for this session list.'};

SL_start         = cfg_entry; 
SL_start.name    = 'Start time(s) for this session list or for each session, in seconds';
SL_start.tag     = 'SL_start';
SL_start.strtype = 's';
SL_start.num     = [0 Inf];
SL_start.val{1}  = '';
SL_start.help    = {'Enter nothing to start at beginning of file, enter 1 number in seconds'
    'to apply to the whole list of session, or specify an array of numbers, on for each session.'};

SL_duration         = cfg_entry; 
SL_duration.name    = 'Duration(s) for this session list or for each session, in seconds';
SL_duration.tag     = 'SL_duration';
SL_duration.strtype = 's';
SL_duration.num     = [0 Inf];
SL_duration.val{1}  = '';
SL_duration.help    = {'Enter nothing to use all available points between specified start and end times, '
    'enter 1 number in seconds to apply to the whole list of session, '
    'or specify an array of numbers, on for each session.'
    'This specifies the duration of the data to use for averaging from the specified start point.'};

SL_end         = cfg_entry; 
SL_end.name    = 'End time(s) for this session list or for each session, in seconds';
SL_end.tag     = 'SL_end';
SL_end.strtype = 's';
SL_end.num     = [0 Inf];
SL_end.val{1}  = '';
SL_end.help    = {'Enter nothing to specify the end of the file, enter 1 number in seconds'
    'to apply to the whole list of session, or specify an array of numbers, on for each session.'
    'This specifies the data not to use for averaging from the end point of the file.'};

session_list2         = cfg_branch;
session_list2.tag     = 'session_list2';
session_list2.name    = 'Session List';
session_list2.val     = {SL_Label SL_List};
session_list2.help    = {'Specify name for list and its content. The list can be a single session if no further averaging is required.'}';

session_list2EXT         = cfg_branch;
session_list2EXT.tag     = 'session_list2EXT';
session_list2EXT.name    = 'Session List (extended)';
session_list2EXT.val     = {SL_Label SL_List SL_start SL_duration SL_end};
session_list2EXT.help    = {'Specify name for list and its content. '
    'The list can be a single session if no further averaging is required.'}';

session_list         = cfg_repeat;
session_list.tag     = 'session_list';
session_list.name    = 'List of sessions';
session_list.help    = {'Create a new item for each list of sessions to be averaged.'}';
session_list.values  = {session_list2};
session_list.num     = [1 Inf];

session_listEXT         = cfg_repeat;
session_listEXT.tag     = 'session_listEXT';
session_listEXT.name    = 'List of sessions (extended)';
session_listEXT.help    = {'Create a new item for each list of sessions to be averaged.'}';
session_listEXT.values  = {session_list2EXT};
session_listEXT.num     = [1 Inf];

average_each_session         = cfg_branch;
average_each_session.tag     = 'average_each_session';
average_each_session.name    = 'Average each session';
average_each_session.val     = {};
average_each_session.help    = {'Simplest option, each session will be averaged separately.'}';

combine_sessions         = cfg_branch;
combine_sessions.tag     = 'combine_sessions';
combine_sessions.name    = 'Average by combining sessions';
combine_sessions.val     = {session_list};
combine_sessions.help    = {'Specify how each session should be treated.'
    'A list must be given, each item of the list being itself a list of sessions to be averaged.'}';

combine_sessionsEXT         = cfg_branch;
combine_sessionsEXT.tag     = 'combine_sessionsEXT';
combine_sessionsEXT.name    = 'Average by combining sessions (extended version)';
combine_sessionsEXT.val     = {session_listEXT};
combine_sessionsEXT.help    = {'Specify how each session should be treated.'
    'A list must be given, each item of the list being itself a list of sessions to be averaged.'
    'In this extended version, one can specify the start and end of each session, and duration too.'
    'Please check within the code that the options work as you intend to.'}';

average_all_data         = cfg_choice;
average_all_data.tag     = 'average_all_data';
average_all_data.name    = 'Average all data';
average_all_data.values  = {average_each_session combine_sessions combine_sessionsEXT};
average_all_data.val     = {average_each_session};
average_all_data.help    = {'Average all data: for each channel this option '
    'will average all time points -- designed for use with the no_baseline option, '
    'and choosing option 3 in module normalize_baseline for the normalization choice'
    '(median of 1st session applied to all sessions).'}';


averaging_choice        = cfg_choice;
averaging_choice.name   = 'Choose averaging method';
averaging_choice.tag    = 'averaging_choice';
averaging_choice.values = {block_averaging average_all_data seizure_evolution};
averaging_choice.val    = {block_averaging};
averaging_choice.help   = {'Choose averaging method.'}';

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vasomotion
%%%%%%%%%%%%%%%%%%%%%%%%%%%

select_chromophore         = cfg_menu;
select_chromophore.tag     = 'select_chromophore';
select_chromophore.name    = 'Select chromophore';
select_chromophore.help    = {'Select chromophore'};
select_chromophore.labels = {'HbT' 'HbR' 'HbO' 'HbO&HbR'};
select_chromophore.values  = {  1   2   3   4};
select_chromophore.val = {1};

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

% ---------------------------------------------------------------------
% sessions Sessions
% ---------------------------------------------------------------------
sessions         = cfg_entry;
sessions.name    = 'Sessions';
sessions.tag     = 'sessions';
sessions.strtype = 'r';
sessions.num     = [0 Inf];
sessions.val     = {''};
sessions.help    = {'Enter the numbers of the sessions to include. '
    'If not specified, all sessions included in the NIRS matrix will be included.'}';

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
time_res.name = 'Time resolution (for downsampling)';
time_res.val = {1};
time_res.strtype = 'r';
time_res.num     = [1 1];
time_res.help    = {'Time resolution for onsets will be given by NIRS'
    'sampling rate divided by this factor  - value is 10 in NIRS_SPM.'}';

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

GenerateHbT      = cfg_menu;
GenerateHbT.tag  = 'GenerateHbT';
GenerateHbT.name = 'Generate HbT';
GenerateHbT.labels = {'Yes','No'};
GenerateHbT.values = {1,0};
GenerateHbT.def = @(val)nirs_get_defaults(...
    'model_specify.wls_bglm_specify.GenerateHbT', val{:});
GenerateHbT.help = {'Generate HbT.'};

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
%hpf_dct_cutoff.def    = @(val)nirs_get_defaults('hpf_dct_cutoff', val{:});
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
nirs_lpf.val       = {lpf_gauss};
nirs_lpf.help      = {'Choose low-pass filter.'};

AvgFilters         = cfg_branch;
AvgFilters.tag     = 'AvgFilters';
AvgFilters.name    = 'AvgFilters';
AvgFilters.val     = {nirs_hpf nirs_lpf};
AvgFilters.help    = {'Specify options for filters.'};

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
hpf_butter_order.help    = {'Enter order of Butterworth HPF (preferred value = 3).'};

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

NumChConfounds         = cfg_entry;
NumChConfounds.name    = 'Maximum Number of Confounds';
NumChConfounds.tag     = 'NumChConfounds';
NumChConfounds.strtype = 'r';
NumChConfounds.num     = [1 1];
NumChConfounds.val     = {1};
NumChConfounds.help    = {'Enter maximum number of NIRS channels to be included as physiological confounds.'};

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

hpf_butter      = cfg_choice;
hpf_butter.tag  = 'hpf_butter';
hpf_butter.name = 'Additional High Pass Filter';
hpf_butter.values = {hpf_butter_On hpf_butter_Off remove_linear GLM_remove_linear SPM_cosine_filter};
hpf_butter.val = {hpf_butter_On};
hpf_butter.help = {'Additional High Pass Filter'}';

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Baseline method
%%%%%%%%%%%%%%%%%%%%%%%%%%%

baseline_offset         = cfg_entry;
baseline_offset.name    = 'Baseline offset in seconds';
baseline_offset.tag     = 'baseline_offset';
baseline_offset.strtype = 'r';
baseline_offset.num     = [1 1];
baseline_offset.val     = {0};
baseline_offset.help    = {'Enter baseline offset in seconds -- typically 0 seconds.'
    'A positive number corresponds to having the end of the baseline window before the stimulus onset.'
    'This option is useful in epilepsy when the hemodynamic response may have'
    'started before the electrophysiological or clinical response.'}';

baseline_duration         = cfg_entry;
baseline_duration.name    = 'Baseline duration in seconds';
baseline_duration.tag     = 'baseline_duration';
baseline_duration.strtype = 'r';
baseline_duration.num     = [1 1];
baseline_duration.val     = {2};
baseline_duration.help    = {'Enter baseline duration.'
    'For an event related protocol, this could be 1 or 2 seconds. '
    'For a block protocol, this will typically be the length of the previous control block.'}';

baseline_block_averaging         = cfg_branch;
baseline_block_averaging.tag     = 'baseline_block_averaging';
baseline_block_averaging.name    = 'Baseline block averaging';
baseline_block_averaging.val     = {baseline_offset baseline_duration};
baseline_block_averaging.help    = {'Baseline block averaging.'
    'This works by averaging the data before each onset over the specified window'}';

baseline_session         = cfg_entry;
baseline_session.name    = 'Baseline session number';
baseline_session.tag     = 'baseline_session';
baseline_session.strtype = 'r';
baseline_session.num     = [1 1];
baseline_session.val     = {2};
baseline_session.help    = {'Enter number of session to be used as baseline.'}';

baseline_block_whole_session         = cfg_branch;
baseline_block_whole_session.tag     = 'baseline_block_whole_session';
baseline_block_whole_session.name    = 'Baseline block averaging using whole session';
baseline_block_whole_session.val     = {baseline_session baseline_offset baseline_duration};
baseline_block_whole_session.help    = {'In this variant, the specified session is used'
    'to compute the baseline for all the other sessions.'
    'Baseline block averaging.'
    'This works by averaging the data before each onset over the specified window'
    'Note that the average will NOT be performed for the session on which the baseline was chosen, contrary to the next option.'}';

baseline_start         = cfg_entry;
baseline_start.name    = 'Baseline start time in seconds';
baseline_start.tag     = 'baseline_start';
baseline_start.strtype = 'r';
baseline_start.num     = [1 1];
baseline_start.val     = {10};
baseline_start.help    = {'Enter start time of baseline.'
    'This will typically be the start of the recording. '
    'However it may be best to remove the first few seconds of the recording '
    'due to edge effects from the high pass filter.'}';

unique_baseline         = cfg_branch;
unique_baseline.tag     = 'unique_baseline';
unique_baseline.name    = 'Unique baseline from chosen session';
unique_baseline.val     = {baseline_session baseline_start baseline_duration};
unique_baseline.help    = {'In this variant, the specified session is used'
    'to compute the baseline for all the other sessions.'
    'This works by subtracting from all onsets the same baseline, typically computed at the beginning of the session.'
    'Note that the average will be performed for the session on which the baseline was chosen, contrary to the previous option.'}';

no_baseline         = cfg_branch;
no_baseline.tag     = 'no_baseline';
no_baseline.name    = 'No baseline';
no_baseline.val     = {};
no_baseline.help    = {'No baseline.'}';

sessions_to_average         = cfg_entry;
sessions_to_average.name    = 'Sessions to average';
sessions_to_average.tag     = 'sessions_to_average';
sessions_to_average.strtype = 'r';
sessions_to_average.num     = [1 Inf];
sessions_to_average.val     = {1};
sessions_to_average.help    = {'Enter single session or list of sessions to average,'
    'this average then to be used as a baseline value.'}';

baseline_average_sessions         = cfg_branch;
baseline_average_sessions.tag     = 'baseline_average_sessions';
baseline_average_sessions.name    = 'For baseline, take average of specified sessions';
baseline_average_sessions.val     = {sessions_to_average SL_start SL_duration SL_end};
baseline_average_sessions.help    = {'For baseline, take average of specified sessions.'}';

baseline_choice        = cfg_choice;
baseline_choice.name   = 'Choose baseline method';
baseline_choice.tag    = 'baseline_choice';
baseline_choice.values = {baseline_block_averaging baseline_block_whole_session ...
    unique_baseline no_baseline baseline_average_sessions};
baseline_choice.val    = {unique_baseline};
baseline_choice.help   = {'Choose baseline method.'}';

% Executable Branch
liom_intrasubject_average      = cfg_exbranch;
liom_intrasubject_average.name = 'LIOM Intrasubject Average';
liom_intrasubject_average.tag  = 'liom_intrasubject_average';
liom_intrasubject_average.val  = {NIRSmat redo1 NIRSmatCopyChoice sessions subj units time_res ...
    GLM_include_cardiac GLM_include_Mayer vasomotion_choice NIRSchannelsConfound GenerateHbT ...
    channel_pca hpf_butter AvgFilters averaging_choice baseline_choice};
liom_intrasubject_average.prog = @nirs_run_liom_intrasubject_average;
liom_intrasubject_average.vout = @nirs_cfg_vout_liom_intrasubject_average;
liom_intrasubject_average.help = {'Specify LIOM General Linear Model.'};

function vout = nirs_cfg_vout_liom_intrasubject_average(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});