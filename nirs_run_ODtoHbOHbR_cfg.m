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
    'as a vector: [PVF(lambda_1) ... PVF(lambda_n)].',...
    'PVF = PPF/DPF (see Strangman 2003, NI 18).'};

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
DPFval.help    = {'Enter the differential pathlength factor values for each wavelength ',...
    'as a vector: [DPF(lambda_1) ... DPF(lambda_n)].',...
    'Default values are for 690 & 830 nm from Duncan 1995.',...
    'PVF = PPF/DPF (see Strangman 2003, NI 18).'};

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

nirs_filling_jumps = nirs_dfg_nirs_filling_jumps;

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