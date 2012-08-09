function ODtoHbOHbR = nirs_run_ODtoHbOHbR_cfg
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'Conc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration for converting Optical Densities to HbO/HbR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PVFsim         = cfg_files;
PVFsim.name    = 'PVF from file';
PVFsim.tag     = 'PVFsim';
PVFsim.filter  = 'mat';
PVFsim.ufilter = '.mat';
PVFsim.num     = [1 1];
PVFsim.help    = {['Select the PVF.mat file created by the nirs_calculatePVE module.']};

PVFval         = cfg_entry;
PVFval.name    = 'Enter values';
PVFval.tag     = 'PVFval';
PVFval.strtype = 'r';
PVFval.num     = [1 2];
PVFval.def  = @(val)nirs_get_defaults('preprocessNIRS.ODtoHbOHbR.PVF', val{:});
PVFval.help    = {'Enter the partial volume factor values for each wavelength ',...
    'as a vector: [PVF(lambda_1) ... PVF(lambda_n)].'};

PVF           = cfg_choice;
PVF.name      = 'Partial Volume Factors';
PVF.tag       = 'PVF';
PVF.values    = {PVFval PVFsim};
PVF.val       = {PVFval};
PVF.help      = {'Either enter a value for the partial volume factor to apply ', ...
    'or use the value computed by Monte-Carlo simulations and stored in the NIRS matrix.'};


DPFsim         = cfg_files;
DPFsim.name    = 'DPF from file';
DPFsim.tag     = 'DPFsim';
DPFsim.filter  = 'mat';
DPFsim.ufilter = '.mat';
DPFsim.num     = [1 Inf];
DPFsim.help    = {['Select the PDPF.mat file(s) created by the nirs_calculatePVE module.',...
    'Select one file or one for each subject in the order corresponding to that in the NIRS matrix.']};

DPFval         = cfg_entry;
DPFval.name    = 'Enter values';
DPFval.tag     = 'DPFval';
DPFval.strtype = 'r';
DPFval.num     = [1 2];
DPFval.def  = @(val)nirs_get_defaults('preprocessNIRS.ODtoHbOHbR.DPF', val{:});
DPFval.help    = {'Enter the partial volume factor values for each wavelength ',...
    'as a vector: [DPF(lambda_1) ... DPF(lambda_n)].'};

DPFlit      = cfg_branch;
DPFlit.name = 'Literature value';
DPFlit.tag  = 'DPFlit';
DPFlit.help = {'Literature values (from Duncan, 1996) will be used.'};

DPF           = cfg_choice;
DPF.name      = 'Differential pathlength factors';
DPF.tag       = 'DPF';
DPF.values    = {DPFval DPFsim DPFlit};
DPF.val       = {DPFlit};
DPF.help      = {'Either enter a value for the differential pathlenth factor to apply ', ...
    'or use the value computed by Monte-Carlo simulations and stored in the NIRS matrix.' ...
    ' Literature values can also be used (default).'};




% ---------------------------------------------------------------------
% lpf Low-pass filter
% ---------------------------------------------------------------------
fwhm2      = cfg_entry;
fwhm2.tag  = 'fwhm2';
fwhm2.name = 'FWHM in seconds';
fwhm2.val = {1.5};
fwhm2.strtype = 'r';
fwhm2.num     = [1 1];
fwhm2.help    = {'FWHM in seconds.'};

downsamplingFactor      = cfg_entry;
downsamplingFactor.tag  = 'downsamplingFactor';
downsamplingFactor.name = 'Downsampling Factor';
downsamplingFactor.val = {1}; %10
downsamplingFactor.strtype = 'r';
downsamplingFactor.num     = [1 1];
downsamplingFactor.help    = {'Specify downsampling factor.'};

downsampleWhen         = cfg_menu;
downsampleWhen.tag     = 'downsampleWhen';
downsampleWhen.name    = 'Apply downsampling on';
downsampleWhen.help    = {'Choose at what step to apply low-pass filtering and downsampling.'};
downsampleWhen.labels = {
    'Raw Optical Densities'
    'Log of Optical Densities'
    'Concentrations'}';
downsampleWhen.values = {1 2 3};
downsampleWhen.val    = {1};

lpf_gauss2         = cfg_branch;
lpf_gauss2.tag     = 'lpf_gauss2';
lpf_gauss2.name    = 'Gaussian Filter';
lpf_gauss2.val     = {fwhm2 downsamplingFactor downsampleWhen};
lpf_gauss2.help    = {'Specify properties of Gaussian filter'};

lpf_none         = cfg_branch;
lpf_none.tag     = 'lpf_none';
lpf_none.name    = 'No low pass filter';
lpf_none.help    = {'No low pass filter.'};

nirs_lpf2           = cfg_choice;
nirs_lpf2.name      = 'Low-pass filter';
nirs_lpf2.tag       = 'nirs_lpf2';
nirs_lpf2.values    = {lpf_none
    lpf_gauss2};
nirs_lpf2.val       = {lpf_none};
nirs_lpf2.help      = {'Choose low-pass filter.'};

%**************************************************************************
%Filling jumps for nirs data
%Ke Peng
%**************************************************************************
hpf_butter_freq         = cfg_entry;
hpf_butter_freq.name    = 'Cutoff frequency for HPF';
hpf_butter_freq.tag     = 'hpf_butter_freq';
hpf_butter_freq.strtype = 'r';
hpf_butter_freq.num     = [1 1];
%hpf_butter_freq.def     = @(val)nirs_get_defaults(...
    %'model_specify.wls_bglm_specify.hpf_butter.hpf_butter_On.hpf_butter_freq', val{:});
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
%cardiac_repair.labels = {'Yes','No'};
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
nirs_filling_jumps.values = {nirs_fill_jumps_on nirs_fill_jumps_off};
nirs_filling_jumps.val = {nirs_fill_jumps_off};
nirs_filling_jumps.help = {'Choose whether to fill the jumps in nirs data.'
                       'Detect jumps and interpolate values to fill them in nirs data'}';
                   
%**************************************************************************


% Executable Branch
ODtoHbOHbR      = cfg_exbranch;
ODtoHbOHbR.name = 'Convert OD to HbO/HbR ';
ODtoHbOHbR.tag  = 'ODtoHbOHbR';
ODtoHbOHbR.val  = {NIRSmat redo1 NIRSmatCopyChoice DPF PVF nirs_filling_jumps}; % nirs_lpf2};
ODtoHbOHbR.prog = @nirs_run_ODtoHbOHbR;
ODtoHbOHbR.vout = @nirs_cfg_vout_ODtoHbOHbR;
ODtoHbOHbR.help = {'Convert OD to HbO/HbR.'}';

function vout = nirs_cfg_vout_ODtoHbOHbR(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});