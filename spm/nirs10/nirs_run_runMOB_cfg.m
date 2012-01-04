function runMOB1 = nirs_run_runMOB_cfg


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Module 14 : CRIUGM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

acc_file         = cfg_files;
acc_file.name    = 'Accelerometer file'; % The displayed name
acc_file.tag     = 'acc_file';       %file names
acc_file.filter  = 'csv';
acc_file.num     = [1 Inf];     % Number of inputs required
acc_file.help    = {''}; % help text displayed

subj         = cfg_branch;
subj.tag     = 'subj';
subj.name    = 'Subject';
subj.val     = {acc_file};
subj.help    = {'Subject'};

generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Subjects';
generic.help    = {'Help'};
generic.values  = {subj};
generic.num     = [1 Inf];

runMOB1      = cfg_exbranch;
runMOB1.name = 'Run MOB analysis';
runMOB1.tag  = 'runMOB1';
runMOB1.val  = {generic};
runMOB1.prog = @nirs_run_runMOB;
runMOB1.vout = @nirs_cfg_vout_runMOB;
runMOB1.help = {'.'};

function vout = nirs_cfg_vout_runMOB(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
