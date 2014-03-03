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

chromophore_select           = cfg_menu;
chromophore_select.name      = 'Select chromophore';
chromophore_select.tag       = 'chromophore_select';
chromophore_select.labels    = {'HbO' 'HbR' 'HbT'};
chromophore_select.values    = {1,2,3};
chromophore_select.val       = {1};
chromophore_select.help      = {'Select a chromophore to estimate HRF'}';

session_select      = cfg_entry;
session_select.tag  = 'session_select';
session_select.name = 'Session';
session_select.num  = [1 1];
session_select.val  = {1};
session_select.help = {'Indicate session number.'};

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
EvtInterest.name = 'Specify the index of event of interest';
EvtInterest.num  = [1 1];
EvtInterest.val  = {0};
EvtInterest.help = {'Specify the number of event of interest. Optional if name if specified below.'};

EvtInterest_name      = cfg_entry;
EvtInterest_name.tag  = 'EvtInterest_name';
EvtInterest_name.name = 'Specify the name of event of interest (optional)';
EvtInterest_name.strtype = 's';
EvtInterest_name.val{1}    = '';
EvtInterest_name.num     = [1 Inf]; 
EvtInterest_name.help = {'Specify the name of event of interest. Optional if index if specified above.'};

EvtNonInterest_sz      = cfg_entry;
EvtNonInterest_sz.tag  = 'EvtNonInterest_sz';
EvtNonInterest_sz.name = 'Specify the index of sz-like events (optional)';
EvtNonInterest_sz.num  = [1 Inf];
EvtNonInterest_sz.val  = {0};
EvtNonInterest_sz.help = {'Specify the index number of all sz-like events. Input 0 if no sz-like events are presented.'};

EvtNonInterest_spk      = cfg_entry;
EvtNonInterest_spk.tag  = 'EvtNonInterest_spk';
EvtNonInterest_spk.name = 'Specify the index of other spk-like events (optional)';
EvtNonInterest_spk.num  = [1 Inf];
EvtNonInterest_spk.val  = {0};
EvtNonInterest_spk.help = {'Specify the index number of other spk-like events. Input 0 if no sz-like events are presented.'};

avoidance_nonintere_only         = cfg_branch;
avoidance_nonintere_only.tag     = 'avoidance_nonintere_only';
avoidance_nonintere_only.name    = 'Only enable avoidance of events of Non-interest';
avoidance_nonintere_only.val     = {EvtNonInterest_sz EvtNonInterest_spk};
avoidance_nonintere_only.help    = {'To enable the avoidance: only avoid the time-period of events of non-interest'};

avoidance_intern_only         = cfg_branch;
avoidance_intern_only.tag     = 'avoidance_intern_only';
avoidance_intern_only.name    = 'Only enable avoidance of interaction of events';
avoidance_intern_only.val     = {};
avoidance_intern_only.help    = {'To enable the avoidance: only avoid the time-period of other events of interest (internal)'};

avoidance_all         = cfg_branch;
avoidance_all.tag     = 'avoidance_all';
avoidance_all.name    = 'Enable all avoidances';
avoidance_all.val     = {EvtNonInterest_sz EvtNonInterest_spk};
avoidance_all.help    = {'To enable all avoidances'};

avoidance_off         = cfg_branch;
avoidance_off.tag     = 'avoidance_off';
avoidance_off.name    = 'Disable avoidance';
avoidance_off.val     = {};
avoidance_off.help    = {'To disable the avoidance'};

TCavoidance        = cfg_choice;
TCavoidance.name   = 'Time course avoidance';
TCavoidance.tag    = 'TCavoidance';
TCavoidance.values = {avoidance_nonintere_only avoidance_intern_only avoidance_all avoidance_off};
TCavoidance.val    = {avoidance_off};
TCavoidance.help   = {'Choose whether to enable time course avoidance.'
    'If enabled, only those events of interest whose time courses are sufficiently apart from those events of non-interest will be included in the analysis'}';

HRFsubtraction           = cfg_menu;
HRFsubtraction.name      = 'Perform subtraction';
HRFsubtraction.tag       = 'HRFsubtraction';
HRFsubtraction.labels    = {'No' 'Yes'};
HRFsubtraction.values    = {0,1};
HRFsubtraction.val       = {0};
HRFsubtraction.help      = {'Select whether subtraction of the averaged HRF is demanded.'}';

simple_averaging         = cfg_branch;
simple_averaging.tag     = 'simple_averaging';
simple_averaging.name    = 'Simple averaging';
simple_averaging.val     = {HRFstart HRFend EvtInterest EvtInterest_name TCavoidance HRFsubtraction};
simple_averaging.help    = {'To use the simple averaging method'};

poly_order      = cfg_entry;
poly_order.tag  = 'poly_order';
poly_order.name = 'Specify the polynomial degree involved in deconvolution';
poly_order.num  = [1 1];
poly_order.val  = {1};
poly_order.help = {'Specify the polynomial degree involved in deconvolution: z(t) = conv(f(t),h(t)) + trend + ipsilon.'};

ols         = cfg_branch;
ols.tag     = 'ols';
ols.name    = 'Ordinary Least Square (OLS)';
ols.val     = {};
ols.help    = {'To use OLS estimators'};

k_val           = cfg_menu;
k_val.name      = 'k value';
k_val.tag       = 'k_val';
k_val.labels    = {'Hoerl single iterative' 'iterative estimation towards convergence'};
k_val.values    = {1,2};
k_val.val       = {1};
k_val.help      = {'Select the method to compute the ridge regression coefficient k.'}';

ridge         = cfg_branch;
ridge.tag     = 'ridge';
ridge.name    = 'Ridge Regression';
ridge.val     = {k_val};
ridge.help    = {'To use ridge regression estimators'};

deriv_kernel         = cfg_choice;
deriv_kernel.tag     = 'deriv_kernel';
deriv_kernel.name    = 'Derivative kernel';
deriv_kernel.values  = {ols ridge};
deriv_kernel.val     = {ols};
deriv_kernel.help    = {'Derivative method: Ordinary Least Square / Ridge regression'};

deconvolution         = cfg_branch;
deconvolution.tag     = 'deconvolution';
deconvolution.name    = 'Deconvolution';
deconvolution.val     = {HRFstart HRFend EvtInterest EvtInterest_name poly_order deriv_kernel};
deconvolution.help    = {'To use the simple averaging method'};

HRFmethod        = cfg_choice;
HRFmethod.name   = 'Method to obtain estimated HRF';
HRFmethod.tag    = 'HRFmethod';
HRFmethod.values = {simple_averaging deconvolution};
HRFmethod.val    = {simple_averaging};
HRFmethod.help   = {'Choose a method to estimate HRF.'}';

channel_no      = cfg_entry;
channel_no.tag  = 'channel_no';
channel_no.name = 'Channel number for output';
channel_no.num  = [1 1];
channel_no.val  = {1};
channel_no.help = {'Specify the desired channel number for output.'};

sample_interval      = cfg_entry;
sample_interval.tag  = 'sample_interval';
sample_interval.name = 'Sample interval (in seconds)';
sample_interval.num  = [1 1];
sample_interval.val  = {1};
sample_interval.help = {'Indicate Sample interval for the estimated HRF'};

channel_view         = cfg_branch;
channel_view.tag     = 'channel_view';
channel_view.name    = 'Channel hemodynamic response variation';
channel_view.val     = {channel_no sample_interval};
channel_view.help    = {'Output the HRF of a certain channel'};

AllowExtrapolation           = cfg_menu;
AllowExtrapolation.name      = 'Allow Extrapolation';
AllowExtrapolation.tag       = 'AllowExtrapolation';
AllowExtrapolation.labels    = {'Yes' 'No'};
AllowExtrapolation.values    = {1,0};
AllowExtrapolation.val       = {0};
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

threshold_value      = cfg_entry;
threshold_value.tag  = 'threshold_value';
threshold_value.name = 'Specify the threshold value';
threshold_value.num  = [1 1];
threshold_value.val  = {2.5};
threshold_value.help = {'Specify the threshold value. Must greater than 0.'};

intepolated_view         = cfg_branch;
intepolated_view.tag     = 'intepolated_view';
intepolated_view.name    = 'Interpolated hemodyanmic response variation';
intepolated_view.val     = {AllowExtrapolation no_interpolation brain_view threshold_value sample_interval};
intepolated_view.help    = {'To use the simple averaging method'};

Outputformat        = cfg_choice;
Outputformat.name   = 'Output format';
Outputformat.tag    = 'Outputformat';
Outputformat.values = {channel_view intepolated_view};
Outputformat.val    = {intepolated_view};
Outputformat.help   = {'Choose output format: interpolated view level / channel level'}';

% Executable Branch
liom_tpHRF      = cfg_exbranch;
liom_tpHRF.name = 'Liom time point HRF';
liom_tpHRF.tag  = 'liom_tpHRF';
liom_tpHRF.val  = {NIRSmat redo1 NIRSmatCopyChoice ...
   chromophore_select session_select HRFmethod Outputformat}; % factorial_design};
liom_tpHRF.prog = @nirs_run_liom_tpHRF;
liom_tpHRF.vout = @nirs_cfg_vout_tpHRF;
liom_tpHRF.help = {'Liom time point HRF estimation.'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vout = nirs_cfg_vout_tpHRF(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
