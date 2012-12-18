function nirs_filling_jumps = nirs_dfg_nirs_filling_jumps
%**************************************************************************
%Filling jumps for nirs data
%Ke Peng
%**************************************************************************
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

HPF_enable_on         = cfg_branch;
HPF_enable_on.tag     = 'HPF_enable_on';
HPF_enable_on.name    = 'Butterworth HPF on';
HPF_enable_on.val     = {hpf_butter_freq hpf_butter_order};
HPF_enable_on.help    = {'Please specify the parameters for Butterworth HPF.'};

HPF_enable_off         = cfg_branch;
HPF_enable_off.tag     = 'HPF_enable_off';
HPF_enable_off.name    = 'Butterworth HPF off';
HPF_enable_off.val     = {};
HPF_enable_off.help    = {'Butterworth HPF has been turned off.'};

HPF_enable      = cfg_choice;
HPF_enable.tag  = 'HPF_enable';
HPF_enable.name = 'Enable Butterworth HPF before filling';
HPF_enable.values = {HPF_enable_on HPF_enable_off};
HPF_enable.val = {HPF_enable_on};
HPF_enable.help = {'Choose whether to enable Butterworth HPF before filling the jumps in nirs data.'};

num_standard_deviation         = cfg_entry;
num_standard_deviation.name    = 'Number of standard deviations:';
num_standard_deviation.tag     = 'num_standard_deviation';
num_standard_deviation.strtype = 'r';
num_standard_deviation.val     = {4};
num_standard_deviation.num     = [1 1];
num_standard_deviation.help    = {'Enter number of standard deviations.'
                    'Integer value must be entered'}';

num_points         = cfg_entry;
num_points.name    = 'Number of points removed before and after:';
num_points.tag     = 'num_points';
num_points.strtype = 'r';
num_points.val     = {1};
num_points.num     = [1 1];
num_points.help    = {'Enter the times of peroid as number of points removed.'
                      'Real value must be entered'}';

size_gap         = cfg_entry;
size_gap.name    = 'Size to define a gap: (in times of period)';
size_gap.tag     = 'size_gap';
size_gap.strtype = 'r';
size_gap.val     = {10};
size_gap.num     = [1 1];
size_gap.help    = {'Enter the times of peroid as a definion of a gap.'
                    'Real value must be entered'}';

nirs_new_fill_jumps         = cfg_branch;
nirs_new_fill_jumps.tag     = 'nirs_new_fill_jumps';
nirs_new_fill_jumps.name    = 'NIRS filling jumps on - new version';
nirs_new_fill_jumps.val     = {};
nirs_new_fill_jumps.help    = {''};

nirs_fill_jumps_on         = cfg_branch;
nirs_fill_jumps_on.tag     = 'nirs_fill_jumps_on';
nirs_fill_jumps_on.name    = 'NIRS filling jumps on';
nirs_fill_jumps_on.val     = {HPF_enable num_standard_deviation num_points size_gap};
nirs_fill_jumps_on.help    = {'Please specify the parameters for gap filling.'};

nirs_fill_jumps_off         = cfg_branch;
nirs_fill_jumps_off.tag     = 'nirs_fill_jumps_off';
nirs_fill_jumps_off.name    = 'NIRS filling jumps off';
nirs_fill_jumps_off.val     = {};
nirs_fill_jumps_off.help    = {'nirs data gap filling has been turned off.'};

nirs_filling_jumps      = cfg_choice;
nirs_filling_jumps.tag  = 'nirs_filling_jumps';
nirs_filling_jumps.name = 'Filling jumps for NIRS data (optional).';
%cardiac_repair.labels = {'Yes','No'};
nirs_filling_jumps.values = {nirs_new_fill_jumps nirs_fill_jumps_on nirs_fill_jumps_off};
nirs_filling_jumps.val = {nirs_fill_jumps_off};
nirs_filling_jumps.help = {'Choose whether to fill the jumps in nirs data.'
                       'Detect jumps and interpolate values to fill them in nirs data'}';
                   
%**************************************************************************
