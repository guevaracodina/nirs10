function liom_cine = nirs_run_liom_cine_cfg
%This is based on nirs_run_liom_contrast, developed to do cinés
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'Cine');
consess = nirs_spm_get_consess;
display_options = liom_contrast_group_options;

%Select view
view         = cfg_entry;
view.name    = 'View';
view.tag     = 'view';
view.strtype = 'r';
view.num     = [1 Inf];
view.val     = {5};
view.help    = {['Enter view.  ',...
    '1: ventral  ',...
    '2: dorsal  ',...
    '3: right  ',...
    '4: left  ',...
    '5: frontal  ',...
    '6: occipital']}; % help text displayed

contrast_p_value         = cfg_entry;
contrast_p_value.name    = 'Corrected p_value';
contrast_p_value.tag     = 'contrast_p_value';
contrast_p_value.strtype = 'r';
contrast_p_value.num     = [1 1];
contrast_p_value.val     = {0.05};
contrast_p_value.help    = {'Corrected p_value'};

spatial_LPF_radius         = cfg_entry;
spatial_LPF_radius.name    = 'Spatial LPF radius';
spatial_LPF_radius.tag     = 'spatial_LPF_radius';
spatial_LPF_radius.strtype = 'r';
spatial_LPF_radius.num     = [1 1];
spatial_LPF_radius.val     = {3};
spatial_LPF_radius.help    = {'Enter radius of spatial low pass filter in pixels.'
    'One pixel is very approximately 1 mm. FWHM will be twice this radius.'
    'If 0 is entered, this is equivalent to no spatial filtering.'
    'Spatial filtering will be applied linearly even though '
    'stereographic projection is nonlinear.'}';

spatial_LPF_On         = cfg_branch;
spatial_LPF_On.tag     = 'spatial_LPF_On';
spatial_LPF_On.name    = 'Spatial LP filter';
spatial_LPF_On.val     = {spatial_LPF_radius};
spatial_LPF_On.help    = {'Spatial low-pass filter.'};

spatial_LPF_Off         = cfg_branch;
spatial_LPF_Off.tag     = 'spatial_LPF_Off';
spatial_LPF_Off.name    = 'Spatial filter off';
spatial_LPF_Off.val     = {};
spatial_LPF_Off.help    = {'Spatial low pass filter turned off.'};

spatial_LPF      = cfg_choice;
spatial_LPF.tag  = 'spatial_LPF';
spatial_LPF.name = 'Spatial Low Pass Filter';
spatial_LPF.values = {spatial_LPF_On spatial_LPF_Off};
spatial_LPF.val = {spatial_LPF_Off};
spatial_LPF.help = {'Choose whether to include a spatial Low Pass Filter'
    'on the interpolated estimates and their variance before constructing'
    'statistical maps.'}';

TopoData        = cfg_files;
TopoData.tag    = 'TopoData';
TopoData.name   = '(Optional) Select TopoData file';
TopoData.filter = '.mat';
TopoData.num    = [0 1];
TopoData.val{1} = {''};
TopoData.help   = {'This is an option for user to specify a new TopoData,'
    'which is generated during coregistration, instead of the default or current TopoData.'}';

InterpolationKernel        = cfg_files;
InterpolationKernel.tag    = 'InterpolationKernel';
InterpolationKernel.name   = '(Optional) Select Interpolation Kernel file';
InterpolationKernel.filter = '.mat';
InterpolationKernel.num    = [0 Inf];
InterpolationKernel.val{1} = {''};
InterpolationKernel.help   = {'This is an option for user to specify a new Interpolation Kernel,'
    'which was generated during a previous run of this cine module.'
    'A list of files can be entered, one for each view,'
    'in the same order as the list of views to be generated.'}';

AllowExtrapolation           = cfg_menu;
AllowExtrapolation.name      = 'Allow Extrapolation';
AllowExtrapolation.tag       = 'AllowExtrapolation';
AllowExtrapolation.labels    = {'Yes' 'No'};
AllowExtrapolation.values    = {1,0};
AllowExtrapolation.val       = {1};
AllowExtrapolation.help      = {'Allow Extrapolation.'}';

no_interpolation           = cfg_menu;
no_interpolation.name      = 'Remove interpolation';
no_interpolation.tag       = 'no_interpolation';
no_interpolation.labels    = {'Yes' 'No'};
no_interpolation.values    = {1,0};
no_interpolation.val       = {0};
no_interpolation.help      = {'Remove interpolation entirely.'}';

% ---------------------------------------------------------------------
% sessions Sessions
% ---------------------------------------------------------------------
Sessions         = cfg_entry; %Careful, not the same as sessions in user-defined SPM contrasts
Sessions.name    = 'Sessions';
Sessions.tag     = 'Sessions';
Sessions.strtype = 'r';
Sessions.num     = [0 Inf];
Sessions.val     = {''};
Sessions.help    = {'Enter the numbers of the sessions to include. If not specified, all sessions present in the NIRS matrix will be included.'};

OnsetChoice         = cfg_entry; %Careful, not the same as sessions in user-defined SPM contrasts
OnsetChoice.name    = 'Onset choice';
OnsetChoice.tag     = 'OnsetChoice';
OnsetChoice.strtype = 'r';
OnsetChoice.num     = [0 Inf];
OnsetChoice.val     = {1};
OnsetChoice.help    = {'Enter the column numbers of the onsets to display.'};

onset_delay         = cfg_entry;
onset_delay.name    = 'Onset delay in seconds';
onset_delay.tag     = 'onset_delay';
onset_delay.strtype = 'r';
onset_delay.num     = [1 1];
onset_delay.val     = {0};
onset_delay.help    = {'Enter onset delay in seconds -- typically 3 seconds.'};

onset_duration         = cfg_entry;
onset_duration.name    = 'Onset duration in seconds';
onset_duration.tag     = 'onset_duration';
onset_duration.strtype = 'r';
onset_duration.num     = [1 1];
onset_duration.val     = {30};
onset_duration.help    = {'Enter onset duration.'
    'For an event related protocol, this could be 1 or 2 seconds. '
    'For a block protocol, this will typically be the length of the block.'}';

time_resolution         = cfg_entry;
time_resolution.name    = 'Time interval between cine frames';
time_resolution.tag     = 'time_resolution';
time_resolution.strtype = 'r';
time_resolution.num     = [1 Inf];
time_resolution.val     = {2};
time_resolution.help    = {'Enter time interval(s) between cine frames, in seconds.'
    'Enter a list to specify non-evenly spaced intervals'}';

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
baseline_duration.val     = {5};
baseline_duration.help    = {'Enter baseline duration.'}';

onsetInfo         = cfg_branch;
onsetInfo.tag     = 'onsetInfo';
onsetInfo.name    = 'Information about onsets'; 
onsetInfo.val     = {OnsetChoice onset_delay onset_duration time_resolution ...
    baseline_offset baseline_duration};
onsetInfo.help    = {'Information about onsets.'}';

select_chromophore         = cfg_entry;
select_chromophore.tag     = 'select_chromophore';
select_chromophore.name    = 'Select chromophore';
select_chromophore.help    = {'Select chromophores: 1 = HbO, 2 = HbR, 3 = HbT'};
select_chromophore.strtype = 'r';
select_chromophore.num     = [1 Inf];
select_chromophore.val     = {[1 2 3]};

%High Pass Filter
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

hpf_butter      = cfg_choice;
hpf_butter.tag  = 'hpf_butter';
hpf_butter.name = 'High Pass Filter';
hpf_butter.values = {hpf_butter_On hpf_butter_Off remove_linear};
hpf_butter.val = {hpf_butter_On};
hpf_butter.help = {'High Pass Filter'}';

% ---------------------------------------------------------------------
% lpf Low-pass filter
% ---------------------------------------------------------------------
fwhm1      = cfg_entry;
fwhm1.tag  = 'fwhm1';
fwhm1.name = 'FWHM in seconds';
fwhm1.val = {1.5};
fwhm1.strtype = 'r';
fwhm1.num     = [1 1];
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

CineFilters         = cfg_branch;
CineFilters.tag     = 'CineFilters';
CineFilters.name    = 'CineFilters';
CineFilters.val     = {hpf_butter nirs_lpf};
CineFilters.help    = {'Specify options for filters.'};

% Executable Branch
liom_cine      = cfg_exbranch;
liom_cine.name = 'Liom CINE';
liom_cine.tag  = 'liom_cine';
liom_cine.val  = {NIRSmat redo1 NIRSmatCopyChoice ...
    view TopoData InterpolationKernel Sessions ...
    select_chromophore onsetInfo ...
    AllowExtrapolation no_interpolation ...
    spatial_LPF CineFilters contrast_p_value display_options}; 
liom_cine.prog = @nirs_run_liom_cine;
liom_cine.vout = @nirs_cfg_vout_liom_cine;
liom_cine.help = {'Liom CINE.'};

function vout = nirs_cfg_vout_liom_cine(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});