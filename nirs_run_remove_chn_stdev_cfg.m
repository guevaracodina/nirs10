function remove_chn_stdev = nirs_run_remove_chn_stdev_cfg

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
%Module 4: Utilities for NIRS data preprocessing: heart rate, filters, etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4.0 Paces... Heart and Mayer -- Mayer not done yet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

threshold_stdev         = cfg_entry;
threshold_stdev.name    = 'Threshold on standard deviation';
threshold_stdev.tag     = 'threshold_stdev';
threshold_stdev.strtype = 'r';
threshold_stdev.num     = [1 1];
threshold_stdev.def     = @(val)nirs_get_defaults('preprocessNIRS.remove_chn_stdev.threshold_stdev', val{:});
threshold_stdev.help    = {'Enter cutoff as a percentage (relative to median of channel) ',...
    'of the median of the standard deviation calculated in rolling windows'.'};

window_stdev         = cfg_entry;
window_stdev.name    = 'Rolling window length';
window_stdev.tag     = 'window_stdev';
window_stdev.strtype = 'r';
window_stdev.num     = [1 1];
window_stdev.def     = @(val)nirs_get_defaults('preprocessNIRS.remove_chn_stdev.window_stdev', val{:});
window_stdev.help    = {'Enter the length in seconds for the rolling windows.'};

% Executable Branch
remove_chn_stdev      = cfg_exbranch;
remove_chn_stdev.name = 'Remove noisy channels (stdev)';
remove_chn_stdev.tag  = 'remove_chn_stdev';
remove_chn_stdev.val  = {NIRSmat DelPreviousData NewDirCopyNIRS threshold_stdev window_stdev};
remove_chn_stdev.prog = @nirs_run_remove_chn_stdev;
remove_chn_stdev.vout = @nirs_cfg_vout_remove_chn_stdev;
remove_chn_stdev.help = {['Preprocessing step: remove noisy channels ',...
    'on the median of a channelwise rolling standard deviation measure.']};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_remove_chn_stdev(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});