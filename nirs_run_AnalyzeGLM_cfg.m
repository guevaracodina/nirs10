function AnalyzeGLM = nirs_run_AnalyzeGLM_cfg
NIRSmat         = cfg_files; %Select NIRS.mat for this subject
NIRSmat.name    = 'NIRS.mat'; % The displayed name
NIRSmat.tag     = 'NIRSmat';       %file names
NIRSmat.filter  = 'mat';
NIRSmat.ufilter = '^NIRS.mat$';
NIRSmat.num     = [1 Inf];     % Number of inputs required
NIRSmat.help    = {'Select NIRS.mat for the subject(s).'}; % help text displayed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyze GLMs - loop over subjects and jobs, to plot simple t contrasts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ROCLoopJob         = cfg_files;
ROCLoopJob.name    = 'Select job(s) to loop over';
ROCLoopJob.tag     = 'ROCLoopJob';
ROCLoopJob.ufilter = '.mat';
ROCLoopJob.num     = [1 Inf];
ROCLoopJob.help    = {'Select .mat-format previously specified '
    'and saved job(s) to loop over.'}';

% Executable Branch
AnalyzeGLM      = cfg_exbranch;
AnalyzeGLM.name = 'Analyze GLMs';
AnalyzeGLM.tag  = 'AnalyzeGLM';
AnalyzeGLM.val  = {NIRSmat ROCLoopJob};
AnalyzeGLM.prog = @nirs_run_AnalyzeGLM;
AnalyzeGLM.vout = @nirs_cfg_vout_AnalyzeGLM;
AnalyzeGLM.help = {'This module performs a large loop over GLMs'
    'from different jobs and subjects.'}';

function vout = nirs_cfg_vout_AnalyzeGLM(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});

