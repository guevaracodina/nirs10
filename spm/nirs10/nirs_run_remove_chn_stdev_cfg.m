function remove_chn_stdev = nirs_run_remove_chn_stdev_cfg
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'std');

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
remove_chn_stdev.val  = {NIRSmat redo1 NIRSmatCopyChoice threshold_stdev window_stdev};
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