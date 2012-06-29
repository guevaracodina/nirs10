function lpf_choice = nirs_dfg_lpf_choice(lpf_on,freq)
% ---------------------------------------------------------------------
% lpf Low-pass filter
% ---------------------------------------------------------------------
fwhm1      = cfg_entry;
fwhm1.tag  = 'fwhm1';
fwhm1.name = 'FWHM in seconds';
fwhm1.val = {freq};
fwhm1.strtype = 'r';
fwhm1.num     = [1 1];
fwhm1.help    = {'FWHM in seconds.'};

lpf_gauss_On         = cfg_branch;
lpf_gauss_On.tag     = 'lpf_gauss_On';
lpf_gauss_On.name    = 'Gaussian LP filter';
lpf_gauss_On.val     = {fwhm1};
lpf_gauss_On.help    = {'Gaussian low-pass filter '}';
% '(applied forward then backward so that it does not create a time shift)'}';

lpf_Off         = cfg_branch;
lpf_Off.tag     = 'lpf_Off';
lpf_Off.name    = 'LP filter off';
lpf_Off.val     = {};
lpf_Off.help    = {'Low pass filter turned off.'};

lpf_choice      = cfg_choice;
lpf_choice.tag  = 'lpf_choice';
lpf_choice.name = 'Choose Low Pass Filter';
lpf_choice.values = {lpf_gauss_On lpf_Off};
if lpf_on
    lpf_choice.val = {lpf_gauss_On};
else
    lpf_choice.val = {lpf_Off};
end
lpf_choice.help = {'Choose whether to include a Low Pass Filter.'
    'Parameters'}';
