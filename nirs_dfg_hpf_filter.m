function hpf_filter = nirs_dfg_hpf_filter
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

cosine_freq         = cfg_entry; 
cosine_freq.name    = 'Cutoff frequency for cosine HPF';
cosine_freq.tag     = 'cosine_freq';       
cosine_freq.strtype = 'r';
cosine_freq.num     = [1 1];     
cosine_freq.val     = {0.01};
cosine_freq.help    = {'Enter cutoff frequency in Hz for SPM Cosine HPF.'};

cosine_freq_band         = cfg_entry; 
cosine_freq_band.name    = 'Cutoff frequency for cosine HPF';
cosine_freq_band.tag     = 'cosine_freq_band';       
cosine_freq_band.strtype = 'r';
cosine_freq_band.num     = [1 2];     
cosine_freq_band.val     = {[0.05 0.12]};
cosine_freq_band.help    = {'Enter one frequency band (e.g. Mayer waves) to remove by '
    'Specifying 2 increasing values'
    'for the range of frequencies to also remove'}';

band_pass_filter_on         = cfg_branch;
band_pass_filter_on.tag     = 'band_pass_filter_on';
band_pass_filter_on.name    = 'Band pass filter on';
band_pass_filter_on.val     = {cosine_freq_band}; 
band_pass_filter_on.help    = {'Band pass filter turned on.'};

band_pass_filter_off         = cfg_branch;
band_pass_filter_off.tag     = 'band_pass_filter_off';
band_pass_filter_off.name    = 'Band pass filter off';
band_pass_filter_off.val     = {}; 
band_pass_filter_off.help    = {'Band pass filter turned off.'};

band_pass_filter         = cfg_choice;
band_pass_filter.tag     = 'band_pass_filter';
band_pass_filter.name    = 'Band pass filter';
band_pass_filter.values  = {band_pass_filter_on band_pass_filter_off};
band_pass_filter.val     = {band_pass_filter_on}; 
band_pass_filter.help    = {'Band pass filte.'};

SPM_cosine_filter         = cfg_branch;
SPM_cosine_filter.tag     = 'SPM_cosine_filter';
SPM_cosine_filter.name    = 'Remove trend using cosines';
SPM_cosine_filter.val     = {cosine_freq band_pass_filter};
SPM_cosine_filter.help    = {'Remove trends using cosines as done in SPM'}';

hpf_filter_Off         = cfg_branch;
hpf_filter_Off.tag     = 'hpf_filter_Off';
hpf_filter_Off.name    = 'HP filter off';
hpf_filter_Off.val     = {}; 
hpf_filter_Off.help    = {'High pass filter turned off.'};

hpf_filter      = cfg_choice;
hpf_filter.tag  = 'hpf_filter';
hpf_filter.name = 'High Pass Filter';
hpf_filter.values = {SPM_cosine_filter hpf_butter_On hpf_filter_Off};
hpf_filter.val = {SPM_cosine_filter};
hpf_filter.help = {'Choose whether to include a DCT filter, or a Butterworth High Pass Filter.'
        'Parameters are: order (e.g. 2) and frequency (e.g. 0.01 Hz), or no filter'}';