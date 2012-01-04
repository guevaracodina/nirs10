function runOT1 = nirs_run_runOvertraining_cfg
NIRSmat         = cfg_files; %Select NIRS.mat for this subject
NIRSmat.name    = 'NIRS.mat'; % The displayed name
NIRSmat.tag     = 'NIRSmat';       %file names
NIRSmat.filter  = 'mat';
NIRSmat.ufilter = '^NIRS.mat$';
NIRSmat.num     = [1 Inf];     % Number of inputs required
NIRSmat.help    = {'Select NIRS.mat for the subject(s).'}; % help text displayed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overtraining ; Olivier Dupuis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


runOT1      = cfg_exbranch;
runOT1.name = 'Run Overtraining analysis';
runOT1.tag  = 'runOT1';
runOT1.val  = {NIRSmat};
runOT1.prog = @nirs_run_runOvertraining;
runOT1.vout = @nirs_cfg_vout_runOvertraining;
runOT1.help = {'.'};

function vout = nirs_cfg_vout_runOvertraining(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
