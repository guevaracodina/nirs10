function ODtoHbOHbR = nirs_run_ODtoHbOHbR_cfg

NIRSmat         = cfg_files; %Select NIRS.mat for this subject
NIRSmat.name    = 'NIRS.mat'; % The displayed name
NIRSmat.tag     = 'NIRSmat';       %file names
NIRSmat.filter  = 'mat';
NIRSmat.ufilter = '^NIRS.mat$';
NIRSmat.num     = [1 Inf];     % Number of inputs required
NIRSmat.help    = {'Select NIRS.mat for the subject(s).'}; % help text displayed

DelPreviousData      = cfg_menu;
DelPreviousData.tag  = 'DelPreviousData';
DelPreviousData.name = 'Delete Previous data file';
DelPreviousData.labels = {'True','False'};
DelPreviousData.values = {1,0};
DelPreviousData.val  = {0};
DelPreviousData.help = {'Delete the previous data file.'}';

CreateNIRSCopy_false         = cfg_branch;
CreateNIRSCopy_false.tag     = 'CreateNIRSCopy_false';
CreateNIRSCopy_false.name    = 'Do not copy NIRS structure';
CreateNIRSCopy_false.help    = {'Do not copy NIRS structure.'
    'This will write over the previous NIRS.mat'}';

NewNIRSdir         = cfg_entry;
NewNIRSdir.name    = 'Directory for NIRS.mat';
NewNIRSdir.tag     = 'NewNIRSdir';
NewNIRSdir.strtype = 's';
NewNIRSdir.val{1}    = 'NewDir';
NewNIRSdir.num     = [1 Inf];
NewNIRSdir.help    = {'Directory for NIRS.mat.'}';

CreateNIRSCopy         = cfg_branch;
CreateNIRSCopy.tag     = 'CreateNIRSCopy';
CreateNIRSCopy.name    = 'Create new directory and copy NIRS structure';
CreateNIRSCopy.val     = {NewNIRSdir};
CreateNIRSCopy.help    = {'Create new directory and copy NIRS structure there.'}';

%Common to most modules: for creating a new directory and copying NIRS.mat
NewDirCopyNIRS           = cfg_choice;
NewDirCopyNIRS.name      = 'Create new directory and copy NIRS.mat';
NewDirCopyNIRS.tag       = 'NewDirCopyNIRS';
NewDirCopyNIRS.values    = {CreateNIRSCopy_false CreateNIRSCopy};
NewDirCopyNIRS.val       = {CreateNIRSCopy_false};
NewDirCopyNIRS.help      = {'Choose whether to overwrite the NIRS.mat structure'
    'or to create a new directory'
    'and copy the NIRS.mat structure there'}';

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
DPFsim.num     = [1 1];
DPFsim.help    = {['Select the PDPF.mat file created by the nirs_calculatePVE module.']};

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


% Executable Branch
ODtoHbOHbR      = cfg_exbranch;
ODtoHbOHbR.name = 'Convert OD to HbO/HbR ';
ODtoHbOHbR.tag  = 'ODtoHbOHbR';
ODtoHbOHbR.val  = {NIRSmat DelPreviousData NewDirCopyNIRS DPF PVF}; % nirs_lpf2};
ODtoHbOHbR.prog = @nirs_run_ODtoHbOHbR;
ODtoHbOHbR.vout = @nirs_cfg_vout_ODtoHbOHbR;
ODtoHbOHbR.help = {'Convert OD to HbO/HbR.'}';

function vout = nirs_cfg_vout_ODtoHbOHbR(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});