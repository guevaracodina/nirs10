function makesens1 = nirs_run_generate_sensitivity_matrix_cfg
NIRSmat         = cfg_files; %Select NIRS.mat for this subject
NIRSmat.name    = 'NIRS.mat'; % The displayed name
NIRSmat.tag     = 'NIRSmat';       %file names
NIRSmat.filter  = 'mat';
NIRSmat.ufilter = '^NIRS.mat$';
NIRSmat.num     = [1 Inf];     % Number of inputs required
NIRSmat.help    = {'Select NIRS.mat for the subject(s).'}; % help text displayed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration: generate sensitivity matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outMCfiles         = cfg_files;
outMCfiles.name    = 'Select MC output files';
outMCfiles.tag     = 'outMCfiles';
outMCfiles.ufilter = {'.2pt','.mc2'};
outMCfiles.num     = [0 Inf];
outMCfiles.val{1}  = {''};
outMCfiles.help    = {'Select .mc2 or .2pt files for this subject.'};

% Executable Branch
makesens1      = cfg_exbranch;
makesens1.name = 'Sensitivity Matrix';
makesens1.tag  = 'makesens1';
makesens1.val  = {NIRSmat outMCfiles};
makesens1.prog = @nirs_run_generate_sensitivity_matrix;
makesens1.vout = @nirs_cfg_vout_generate_sensitivity_matrix;
makesens1.help = {'Generate sensitivity matrix.'};

function vout = nirs_cfg_vout_generate_sensitivity_matrix(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});