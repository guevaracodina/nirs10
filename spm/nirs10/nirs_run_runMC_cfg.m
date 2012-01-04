function runMC1 = nirs_run_runMC_cfg

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
%Configuration: run MC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MCXconfigFiles         = cfg_files;
MCXconfigFiles.name    = 'Select input files';
MCXconfigFiles.tag     = 'MCXconfigFiles';
MCXconfigFiles.ufilter = '.inp';
MCXconfigFiles.num     = [0 Inf];
MCXconfigFiles.val{1}  = {''};
MCXconfigFiles.help    = {'Select input files (.inp for MCX).'};

%%%%% Options for MC run
MCX_t      = cfg_entry;
MCX_t.tag  = 'MCX_t';
MCX_t.name = 'Thread number';
MCX_t.val = {4800};
MCX_t.strtype = 'r';
MCX_t.num     = [1 1];
MCX_t.help = {'Total number of threads -- see the specifications of your GPU'
    'Examples: NVidia GeForce 570: 4800'
    'GeForce GTX295: 1792?'}';

MCX_T      = cfg_entry;
MCX_T.tag  = 'MCX_T';
MCX_T.name = 'Thread number per block';
MCX_T.val = {480};
MCX_T.strtype = 'r';
MCX_T.num     = [1 1];
MCX_T.help = {'Blocksize -- see the specifications of your GPU'
    'Examples: NVidia GeForce 570: 480'
    'GeForce GTX295: ?'}';

MCX_r      = cfg_entry;
MCX_r.tag  = 'MCX_r';
MCX_r.name = 'Number of repetitions';
MCX_r.val = {1};
MCX_r.strtype = 'r';
MCX_r.num     = [1 1];
MCX_r.help = {'Number of repetitions: number of times that the simulation'
    'will be repeated, with different random seeds, to increase the total'
    'number of photons'}';

MCX_g      = cfg_entry;
MCX_g.tag  = 'MCX_g';
MCX_g.name = 'Number of gates';
MCX_g.val = {1};
MCX_g.strtype = 'r';
MCX_g.num     = [1 1];
MCX_g.help = {'Number of gates: if larger than the number of specified gates'
    'in the config files, then only the number of gates specified in the config'
    'files will be run. Otherwise, this allows the user to run with fewer gates.'}';

MCX_l      = cfg_menu;
MCX_l.tag  = 'MCX_l';
MCX_l.name = 'Write log file';
MCX_l.labels = {'Yes', 'No'};
MCX_l.values = {1,0};
MCX_l.val  = {1};
MCX_l.help = {'Write log file.'}';

MCX1         = cfg_branch;
MCX1.tag     = 'MCX1';
MCX1.name    = 'Monte Carlo Extreme';
MCX1.val     = {MCXconfigFiles MCX_t MCX_T MCX_r MCX_g MCX_l}; % MCXconfig};
MCX1.help    = {'Run Monte Carlo Extreme simulation'};

tMCimg_configFiles         = cfg_files; %Select
tMCimg_configFiles.name    = 'Select input files';
tMCimg_configFiles.tag     = 'tMCimg_configFiles';
tMCimg_configFiles.ufilter = '.cfg';
tMCimg_configFiles.num     = [0 Inf];
tMCimg_configFiles.val{1}  = {''};
tMCimg_configFiles.help    = {'Select input files (.cfg for tMCimg).'};

tMCimg1         = cfg_branch;
tMCimg1.tag     = 'tMCimg1';
tMCimg1.name    = 'tMCimg Monte Carlo Simulation';
tMCimg1.val     = {tMCimg_configFiles};
tMCimg1.help    = {'Run tMCimg Monte Carlo simulation'};
%
MC_runCUDAchoice        = cfg_choice;
MC_runCUDAchoice.name   = 'Monte Carlo simulation method';
MC_runCUDAchoice.tag    = 'MC_runCUDAchoice';
MC_runCUDAchoice.values = {MCX1,tMCimg1};
MC_runCUDAchoice.val    = {MCX1};
MC_runCUDAchoice.help   = {['Choose method of Monte Carlo simulation. ',...
    'MCX is much faster but requires a CUDA compatible graphics card.']};

MCtestOneChannel      = cfg_menu;
MCtestOneChannel.tag  = 'MCtestOneChannel';
MCtestOneChannel.name = 'Test MC by running only first channel';
MCtestOneChannel.labels = {'Yes','No'};
MCtestOneChannel.values = {1,0};
MCtestOneChannel.val  = {0};
MCtestOneChannel.help = {'To do a quick test, run simulation only on first source and first detector for one wavelength.'}';


% Executable Branch
runMC1      = cfg_exbranch;
runMC1.name = 'Run Monte Carlo simulation';
runMC1.tag  = 'runMC1';
runMC1.val  = {NIRSmat NewDirCopyNIRS MC_runCUDAchoice MCtestOneChannel};
runMC1.prog = @nirs_run_runMC;
runMC1.vout = @nirs_cfg_vout_runMC;
runMC1.help = {'Run Monte Carlo simulation.'};

%make .mc2 or (.his, .2pt) file names available as a dependency
function vout = nirs_cfg_vout_runMC(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});