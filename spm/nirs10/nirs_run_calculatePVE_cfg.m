function calculatePVE1 = nirs_run_calculatePVE_cfg
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'PVE');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4.8 Calculate partial volume effect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_in         = cfg_files;
dir_in.tag     = 'dir_in';
dir_in.name    = 'MonteCarlo output directory';
dir_in.help    = {'Select the MonteCarlo simulation output directory.'};
dir_in.filter = 'dir';
dir_in.val{1} = {''};
dir_in.ufilter = '.*';
dir_in.num     = [0 1];

% historyfiles         = cfg_files;
% historyfiles.name    = 'Monte Carlo history files';
% historyfiles.tag     = 'historyfiles';
% historyfiles.ufilter = {'.his','.mch'};
% historyfiles.num     = [1 Inf];
% historyfiles.help    = {'Select history files for this subject.'};

% Executable Branch
calculatePVE1      = cfg_exbranch;
calculatePVE1.name = 'Calculate Partial Volume Effect';
calculatePVE1.tag  = 'calculatePVE1';
calculatePVE1.val  = {NIRSmat redo1 NIRSmatCopyChoice dir_in};
calculatePVE1.prog = @nirs_run_calculatePVE;
calculatePVE1.vout = @nirs_cfg_vout_calculatePVE;
calculatePVE1.help = {'Calculate Partial Volume Effect'};

function vout = nirs_cfg_vout_calculatePVE(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});