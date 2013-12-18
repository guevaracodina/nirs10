function liom_tpHRF = nirs_run_liom_tpHRF_cfg

[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'tpHRF');

TopoData        = cfg_files;
TopoData.tag    = 'TopoData';
TopoData.name   = '(Optional) Select TopoData file';
TopoData.filter = '.mat';
TopoData.num    = [0 1];
TopoData.val{1} = {''};
TopoData.help   = {'This is an option for user to specify a new TopoData,'
    'which is generated during coregistration, instead of the default or current TopoData.'}';

AllowExtrapolation           = cfg_menu;
AllowExtrapolation.name      = 'Allow Extrapolation';
AllowExtrapolation.tag       = 'AllowExtrapolation';
AllowExtrapolation.labels    = {'Yes' 'No'};
AllowExtrapolation.values    = {1,0};
AllowExtrapolation.val       = {1};
AllowExtrapolation.help      = {'Allow Extrapolation.'}';

no_interpolation           = cfg_menu;
no_interpolation.name      = 'Remove interpolation';
no_interpolation.tag       = 'no_interpolation';
no_interpolation.labels    = {'Yes' 'No'};
no_interpolation.values    = {1,0};
no_interpolation.val       = {0};
no_interpolation.help      = {'Remove interpolation entirely.'}';

brain_view           = cfg_menu;
brain_view.name      = 'Select a brain view';
brain_view.tag       = 'brain_view';
brain_view.labels    = {'Dorsal' 'Right' 'Left' 'Frontal' 'Occipital'};
brain_view.values    = {2,3,4,5,6};
brain_view.val       = {2};
brain_view.help      = {'Select a brain view to create map.'}';

chromophore_select           = cfg_menu;
chromophore_select.name      = 'Select chromophore';
chromophore_select.tag       = 'chromophore_select';
chromophore_select.labels    = {'HbO' 'HbR' 'HbT'};
chromophore_select.values    = {0,1,2};
chromophore_select.val       = {0};
chromophore_select.help      = {'Select a chromophore to estimate HRF'}';

session_select      = cfg_entry;
session_select.tag  = 'session_select';
session_select.name = 'Session';
session_select.num  = [1 1];
session_select.val  = {1};
session_select.help = {'Indicate session number.'};

sample_interval      = cfg_entry;
sample_interval.tag  = 'sample_interval';
sample_interval.name = 'Sample interval (in seconds)';
sample_interval.num  = [1 1];
sample_interval.val  = {1};
sample_interval.help = {'Indicate Sample interval for the estimated HRF'};

HRFstart      = cfg_entry;
HRFstart.tag  = 'HRFstart';
HRFstart.name = 'HRF start time (in seconds)';
HRFstart.num  = [1 1];
HRFstart.val  = {-20};
HRFstart.help = {'Specify the start time of persumed HRF, both positive number and negative number are accepted'};

HRFend      = cfg_entry;
HRFend.tag  = 'HRFend';
HRFend.name = 'HRF end time (in seconds)';
HRFend.num  = [1 1];
HRFend.val  = {30};
HRFend.help = {'Specify the end time of persumed HRF, both positive number and negative number are accepted'};

EvtInterest      = cfg_entry;
EvtInterest.tag  = 'EvtInterest';
EvtInterest.name = 'Specify the event of interest (optional)';
EvtInterest.num  = [1 1];
EvtInterest.val  = {1};
EvtInterest.help = {'Specify the number of event of interest. If non-specified, the first event will be taken.'};

EvtNonInterest      = cfg_entry;
EvtNonInterest.tag  = 'EvtNonInterest';
EvtNonInterest.name = 'Specify the event of Non interest (optional)';
EvtNonInterest.num  = [1 Inf];
EvtNonInterest.val  = {2};
EvtNonInterest.help = {'Specify the number of event of interest. If non-specified, the second event will be taken.'};

avoidance_on         = cfg_branch;
avoidance_on.tag     = 'avoidance_on';
avoidance_on.name    = 'Enable avoidance';
avoidance_on.val     = {EvtNonInterest};
avoidance_on.help    = {'To enable the avoidance'};

avoidance_off         = cfg_branch;
avoidance_off.tag     = 'avoidance_off';
avoidance_off.name    = 'Disable avoidance';
avoidance_off.val     = {};
avoidance_off.help    = {'To disable the avoidance'};

TCavoidance        = cfg_choice;
TCavoidance.name   = 'Time course avoidance';
TCavoidance.tag    = 'TCavoidance';
TCavoidance.values = {avoidance_on avoidance_off};
TCavoidance.val    = {avoidance_on};
TCavoidance.help   = {'Choose whether to enable time course avoidance.'
    'If enabled, only those events of interest whose time courses are sufficiently apart from those events of non-interest will be included in the analysis'}';

simple_averaging         = cfg_branch;
simple_averaging.tag     = 'simple_averaging';
simple_averaging.name    = 'Simple averaging';
simple_averaging.val     = {HRFstart HRFend EvtInterest TCavoidance};
simple_averaging.help    = {'To use the simple averaging method'};

HRFmethod        = cfg_choice;
HRFmethod.name   = 'Method to obtain estimated HRF';
HRFmethod.tag    = 'HRFmethod';
HRFmethod.values = {simple_averaging};
HRFmethod.val    = {simple_averaging};
HRFmethod.help   = {'Choose a method to estimate HRF.'}';

% Executable Branch
liom_tpHRF      = cfg_exbranch;
liom_tpHRF.name = 'Liom time point HRF';
liom_tpHRF.tag  = 'liom_tpHRF';
liom_tpHRF.val  = {NIRSmat redo1 NIRSmatCopyChoice ...
   AllowExtrapolation no_interpolation brain_view chromophore_select session_select sample_interval HRFmethod}; % factorial_design};
liom_tpHRF.prog = @nirs_run_liom_tpHRF;
liom_tpHRF.vout = @nirs_cfg_vout_tpHRF;
liom_tpHRF.help = {'Liom time point HRF estimation.'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vout = nirs_cfg_vout_tpHRF(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
