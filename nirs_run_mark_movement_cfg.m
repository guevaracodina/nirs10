function mark_movement = nirs_run_mark_movement_cfg
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'mvt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Preprocess NIRS: mark movement jumps and shifts as bad and normalize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mvt_ch_thresh         = cfg_entry;
mvt_ch_thresh.name    = 'Movement time cutoff';
mvt_ch_thresh.tag     = 'mvt_ch_thresh';
mvt_ch_thresh.strtype = 'r';
mvt_ch_thresh.num     = [1 1];
mvt_ch_thresh.def     = @(val)nirs_get_defaults('preprocessNIRS.mark_movement.mvt_ch_thresh', val{:});
mvt_ch_thresh.help    = {'Enter the maximal percentage of time allowed '
    'for a good channel - this will be tested in the first session'
    'channels that failed will be excluded for this and following sessions.'}';

mvt_window_length         = cfg_entry;
mvt_window_length.name    = 'Window size to detect movement';
mvt_window_length.tag     = 'mvt_window_length';
mvt_window_length.strtype = 'r';
mvt_window_length.num     = [1 1];
mvt_window_length.def     = @(val)nirs_get_defaults('preprocessNIRS.mark_movement.mvt_window_length', val{:});
mvt_window_length.help    = {['Enter the length of the window in seconds ',...
    'over which to detect movement.']};

mvt_cutoff         = cfg_entry;
mvt_cutoff.name    = 'Movement cutoff';
mvt_cutoff.tag     = 'mvt_cutoff';
mvt_cutoff.strtype = 'r';
mvt_cutoff.num     = [1 1];
mvt_cutoff.def     = @(val)nirs_get_defaults('preprocessNIRS.mark_movement.mvt_cutoff', val{:});
mvt_cutoff.help    = {'Enter the maximal change allowed '
    'as a percentage of the median intensity signal over the window '
    'length, above which data points over the window length are marked as bad.'}';

sum_mvt_threshold         = cfg_entry;
sum_mvt_threshold.name    = 'Threshold on number of channels';
sum_mvt_threshold.tag     = 'sum_mvt_threshold';
sum_mvt_threshold.strtype = 'r';
sum_mvt_threshold.num     = [1 1];
sum_mvt_threshold.def     = @(val)nirs_get_defaults('preprocessNIRS.mark_movement.sum_mvt_threshold', val{:});
sum_mvt_threshold.help    = {'Enter the value of the threshold as '
    'a percentage of the total number of channels, so that '
    'data points which are marked bad for at least that many channels '
    'will be marked bad for all channels.'}';

min_session_duration         = cfg_entry;
min_session_duration.name    = 'Minimum length of subsessions';
min_session_duration.tag     = 'min_session_duration';
min_session_duration.strtype = 'r';
min_session_duration.num     = [1 1];
min_session_duration.def     = @(val)nirs_get_defaults('preprocessNIRS.mark_movement.min_session_duration', val{:});
min_session_duration.help    = {'Enter the minimum length (for example 60 seconds).'
    'of subsessions, i.e. the minimum intervals between '
    'markers of movement.'}';

% Executable Branch
mark_movement      = cfg_exbranch;
mark_movement.name = 'Mark movement';
mark_movement.tag  = 'mark_movement';
mark_movement.val  = {NIRSmat redo1 NIRSmatCopyChoice mvt_ch_thresh...
    mvt_window_length mvt_cutoff sum_mvt_threshold min_session_duration};
mark_movement.prog = @nirs_run_mark_movement;
mark_movement.vout = @nirs_cfg_vout_mark_movement;
mark_movement.help = {['Preprocessing step: mark movement jumps ',...
    'and shifts, segment into intervals.']};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_mark_movement(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});