function ROCtest = nirs_run_ROCtest_cfg

NIRSmat         = cfg_files; %Select NIRS.mat for this subject
NIRSmat.name    = 'NIRS.mat'; % The displayed name
NIRSmat.tag     = 'NIRSmat';       %file names
NIRSmat.filter  = 'mat';
NIRSmat.ufilter = '^NIRS.mat$';
NIRSmat.num     = [1 Inf];     % Number of inputs required
NIRSmat.help    = {'Select NIRS.mat for the subject(s).'}; % help text displayed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROC - Receiver Operating Curve Module - Sensitivity and specificity test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ROCDeleteLarge      = cfg_menu;
ROCDeleteLarge.tag  = 'ROCDeleteLarge';
ROCDeleteLarge.name = 'Delete large files';
ROCDeleteLarge.labels = {'Yes','No'};
ROCDeleteLarge.values = {1,2};
ROCDeleteLarge.val = {1};
ROCDeleteLarge.help = {'Delete large files (.nir) after each estimation.'};

ROCiternum         = cfg_entry;
ROCiternum.name    = 'Number of iterations';
ROCiternum.tag     = 'ROCiternum';
ROCiternum.val     = {10};
ROCiternum.strtype = 'r';
ROCiternum.num     = [1 1];
ROCiternum.help    = {'Number of iterations'};

RunGLMorFigures      = cfg_menu;
RunGLMorFigures.tag  = 'RunGLMorFigures';
RunGLMorFigures.name = 'Run GLMs or generate figures';
RunGLMorFigures.labels = {'GLM','Figures','Both'};
RunGLMorFigures.values = {1,2,3};
RunGLMorFigures.val    = {3};
RunGLMorFigures.help = {'Run GLMs and/or generate figures from'
    'Previously run GLMs.'}';

%ROC options
ROCnumCh         = cfg_entry;
ROCnumCh.name    = 'Total Number of channels';
ROCnumCh.tag     = 'ROCnumCh';
ROCnumCh.val     = {40};
ROCnumCh.strtype = 'r';
ROCnumCh.num     = [1 1];
ROCnumCh.help    = {'Total Number of channels'};

dir_dataSPM         = cfg_entry;
dir_dataSPM.name    = 'Directory to work from';
dir_dataSPM.tag     = 'dir_dataSPM';
dir_dataSPM.val{1}  = 'dataSPM';
dir_dataSPM.strtype = 's';
dir_dataSPM.num     = [1 Inf];
dir_dataSPM.help    = {'Directory to work from'};

Volt2         = cfg_entry;
Volt2.name    = 'Positive or negative t-test for 2nd Volterra';
Volt2.tag     = 'Volt2';
Volt2.val     = {0};
Volt2.strtype = 'r';
Volt2.num     = [1 Inf];
Volt2.help    = {'Positive or negative t-test for 2nd Volterra'
    'Enter a matrix of number of jobs by number of subjects'
    'With entries of 0 if negative t-test and 1 if positive t-test'
    'for the 2nd Volterra'
    'Or enter just 0 (and not a matrix) if all t-tests are negative.'}';

byIter           = cfg_menu;
byIter.name      = 'Give test result by iteration';
byIter.tag       = 'byIter';
byIter.labels    = {'No' 'Yes'};
byIter.values    = {0,1};
byIter.val       = {0};
byIter.help      = {'Usually, No.'}';

compute_OR           = cfg_menu;
compute_OR.name      = 'Compute HbO and HbR separately';
compute_OR.tag       = 'compute_OR';
compute_OR.labels    = {'No' 'Yes'};
compute_OR.values    = {0,1};
compute_OR.val       = {0};
compute_OR.help      = {'Usually, No.'}';

compute_LU           = cfg_menu;
compute_LU.name      = 'Compute lower and upper bounds';
compute_LU.tag       = 'compute_LU';
compute_LU.labels    = {'No' 'Yes'};
compute_LU.values    = {0,1};
compute_LU.val       = {0};
compute_LU.help      = {'Usually, No.'}';

runFtest           = cfg_menu;
runFtest.name      = 'Run F test';
runFtest.tag       = 'runFtest';
runFtest.labels    = {'No' 'Yes'};
runFtest.values    = {0,1};
runFtest.val       = {0};
runFtest.help      = {'Usually, No.'
    'Careful, can be very slow, since...'}';


ROCLoopJob         = cfg_files;
ROCLoopJob.name    = 'Select job(s) to loop over';
ROCLoopJob.tag     = 'ROCLoopJob';
ROCLoopJob.ufilter = '.mat';
ROCLoopJob.num     = [1 Inf];
ROCLoopJob.help    = {'Select .mat-format previously specified '
    'and saved job(s) to loop over.'}';


% Executable Branch
ROCtest      = cfg_exbranch;
ROCtest.name = 'ROC Sensitivity and specificity testing';
ROCtest.tag  = 'ROCtest';
ROCtest.val  = {NIRSmat ROCLoopJob ROCDeleteLarge ROCiternum RunGLMorFigures ...
    ROCnumCh dir_dataSPM Volt2 byIter compute_OR compute_LU runFtest};
ROCtest.prog = @nirs_run_ROCtest;
ROCtest.vout = @nirs_cfg_vout_ROCtest;
ROCtest.help = {'This module performs a large loop over GLMs'
    'specified with different random seeds. '
    'To use it, user need to first specify in the Matlabbatch'
    'front end a sequence of modules to be run, typically starting'
    'with a module that requires a random seed (such as the '
    'AddTestStimuli module). This sequence of modules is referred '
    'to as a job. The code will run that job repetitively '
    'by incrementing the random seed as many times as specified. '
    'One can later use the separate script ROCfigures_script'
    'to generate a variety of figures: ROC plots and boxplots.'
    'This script will require adapting it to specific requirements.'}';

function vout = nirs_cfg_vout_ROCtest(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});