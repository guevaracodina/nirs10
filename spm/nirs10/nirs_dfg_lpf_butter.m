function lpf_butter = nirs_dfg_lpf_butter
lpf_butter_freq         = cfg_entry;
lpf_butter_freq.name    = 'Cutoff frequency for LPF';
lpf_butter_freq.tag     = 'lpf_butter_freq';
lpf_butter_freq.strtype = 'r';
lpf_butter_freq.num     = [1 1];
lpf_butter_freq.val     = {0.667};
lpf_butter_freq.help    = {'Enter cutoff frequency in Hz for Butterworth LPF.'};

lpf_butter_order         = cfg_entry;
lpf_butter_order.name    = 'Order of Butterworth LPF';
lpf_butter_order.tag     = 'lpf_butter_order';
lpf_butter_order.strtype = 'r';
lpf_butter_order.num     = [1 1];
lpf_butter_order.val     = {4};
lpf_butter_order.help    = {'Enter order of Butterworth LPF (preferred value = 2).'};

lpf_butter         = cfg_branch;
lpf_butter.tag     = 'lpf_butter';
lpf_butter.name    = 'Butterworth LP filter';
lpf_butter.val     = {lpf_butter_freq lpf_butter_order};
lpf_butter.help    = {'Butterworth low-pass filter.'};