function makesens1 = nirs_run_generate_sensitivity_matrix_cfg
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'sens');

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
makesens1.val  = {NIRSmat redo1 NIRSmatCopyChoice outMCfiles};
makesens1.prog = @nirs_run_generate_sensitivity_matrix;
makesens1.vout = @nirs_cfg_vout_generate_sensitivity_matrix;
makesens1.help = {'Generate sensitivity matrix.'};

function vout = nirs_cfg_vout_generate_sensitivity_matrix(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});