function coreg_manual1 = nirs_run_coreg_manual_cfg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration for coregistration: coreg MANUAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


NIRSmatSingle       = cfg_files; %Select NIRS.mat for this subject
NIRSmatSingle.name    = 'Select NIRS.mat'; % The displayed name
NIRSmatSingle.tag     = 'NIRSmatSingle';       %file names
NIRSmatSingle.filter  = 'mat';
NIRSmatSingle.ufilter = '^NIRS.mat$';
NIRSmatSingle.num     = [1 1];     % Number of inputs required
NIRSmatSingle.help    = {'Select NIRS.mat for this subject.'};

CoregFromNIRS         = cfg_branch;
CoregFromNIRS.tag     = 'CoregFromNIRS';
CoregFromNIRS.name    = 'Manual coregistration using NIRS.mat';
CoregFromNIRS.val     = {NIRSmatSingle}; % tMCimg_config};
CoregFromNIRS.help    = {'Manual coregistration using previously created NIRS.mat'};

WanatT1         = cfg_files; %Select T1 for this subject
WanatT1.name    = 'Select normalized anatomical image';
WanatT1.tag     = 'WanatT1';       %file names
WanatT1.filter  = 'image';
WanatT1.ufilter = '.*';
WanatT1.num     = [1 1];     % Number of inputs required
WanatT1.help    = {'Select normalized anatomical image for this subject.'};

NormParams         = cfg_files; %Select T1 for this subject
NormParams.name    = 'Select normalized parameters';
NormParams.tag     = 'NormParams';       %file names
%NormParams.filter  = '';
NormParams.ufilter = '_sn.*';
NormParams.num     = [1 1];     % Number of inputs required
NormParams.help    = {'Select normalization parameters for this subject.'};


Coreg_standalone         = cfg_branch;
Coreg_standalone.tag     = 'Coreg_standalone';
Coreg_standalone.name    = 'Standalone coregistration';
Coreg_standalone.val     = {WanatT1 NormParams}; % tMCimg_config};
Coreg_standalone.help    = {'Select normalized anatomical image for optode positioning for this subject'};

Coreg_choice        = cfg_choice;
Coreg_choice.name   = 'Manual coregistration choice';
Coreg_choice.tag    = 'Coreg_choice';
Coreg_choice.values = {CoregFromNIRS,Coreg_standalone};
Coreg_choice.val    = {CoregFromNIRS};
Coreg_choice.help   = {'Choose whether to run from NIRS.mat or as standalone'};

Vsegmented         = cfg_files; %Select anatomical image for this subject
Vsegmented.name    = 'Select anatomical image'; % The displayed name
Vsegmented.tag     = 'Vsegmented';       %file names
Vsegmented.filter  = 'image';
Vsegmented.ufilter = '.*';
Vsegmented.num     = [1 1];     % Number of inputs required
Vsegmented.help    = {'Select native (not normalized) anatomical image for this subject.'};

% Executable Branch
coreg_manual1      = cfg_exbranch;
coreg_manual1.name = 'Manual NIRScoreg';
coreg_manual1.tag  = 'coreg_manual1';
coreg_manual1.val  = {Coreg_choice Vsegmented };
coreg_manual1.prog = @nirs_run_coreg_manual;
coreg_manual1.vout = @nirs_cfg_vout_coreg_manual;
coreg_manual1.help = {'Manual coregistration.'};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_coreg_manual(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});