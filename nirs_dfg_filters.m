function filters = nirs_dfg_filters

channel_pca      = cfg_menu;
channel_pca.tag  = 'channel_pca';
channel_pca.name = 'Spatial Principal Component Removal';
channel_pca.labels = {'Yes','No'};
channel_pca.values = {1,0};
channel_pca.val{1} = 0;
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
% 
% fwhm1      = cfg_entry;
% fwhm1.tag  = 'fwhm1';
% fwhm1.name = 'FWHM in seconds';
% fwhm1.val = {1.5};
% fwhm1.strtype = 'r';
% fwhm1.num     = [1 1];
% fwhm1.help    = {'FWHM in seconds.'};

% lpf_gauss         = cfg_branch;
% lpf_gauss.tag     = 'lpf_gauss';
% lpf_gauss.name    = 'Gaussian Filter';
% lpf_gauss.val     = {fwhm1};
% lpf_gauss.help    = {'Specify properties of Gaussian filter'};

lpf_none         = cfg_branch;
lpf_none.tag     = 'lpf_none';
lpf_none.name    = 'No low pass filter';
lpf_none.help    = {'No low pass filter.'};

% lpf_hrf         = cfg_branch;
% lpf_hrf.tag     = 'lpf_hrf';
% lpf_hrf.name    = 'HRF Filter';
% lpf_hrf.help    = {'HRF filter'};
% 
% lpf           = cfg_choice;
% lpf.name      = 'Low-pass filter';
% lpf.tag       = 'lpf';
% lpf.values    = {lpf_none lpf_gauss lpf_hrf};
% lpf.val       = {lpf_hrf};
% lpf.help      = {'Choose low-pass filter.'}';

lpf_butter = nirs_dfg_lpf_butter;

lpf           = cfg_choice;
lpf.name      = 'Low-pass filter';
lpf.tag       = 'lpf';
lpf.values    = {lpf_none lpf_butter};
lpf.val       = {lpf_butter};
lpf.help      = {'Choose low-pass filter.'}';

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

SPM_cosine_filter         = cfg_branch;
SPM_cosine_filter.tag     = 'SPM_cosine_filter';
SPM_cosine_filter.name    = 'Remove trend using cosines in the GLM';
SPM_cosine_filter.val     = {};
SPM_cosine_filter.help    = {'Remove trends using cosines as done in SPM'}';

hpf      = cfg_choice;
hpf.tag  = 'hpf';
hpf.name = 'Additional High Pass Filter';
hpf.values = {hpf_butter_On hpf_butter_Off remove_linear SPM_cosine_filter};
hpf.val = {hpf_butter_On};
hpf.help = {'Additional High Pass Filter'}';

filters_on      = cfg_branch;
filters_on.tag  = 'filters_on';
filters_on.name = 'Include filters on raw data';
filters_on.val = {channel_pca NumPCAComponents lpf hpf};
filters_on.help = {'Include filters'}';

filters_off      = cfg_branch;
filters_off.tag  = 'filters_off';
filters_off.name = 'Do not include filters on raw data';
filters_off.val = {};
filters_off.help = {'Do not include filters on raw data'}';

filters      = cfg_choice;
filters.tag  = 'filters';
filters.name = 'Include filters directly on raw data';
filters.values = {filters_off filters_on};
filters.val = {filters_off};
filters.help = {'Include filters directly on raw data'}';




