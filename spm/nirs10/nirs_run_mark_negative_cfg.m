function mark_negative = nirs_run_mark_negative_cfg
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'neg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Preprocess NIRS: mark negative points as bad and interpolate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sum_neg_threshold         = cfg_entry;
sum_neg_threshold.name    = 'Threshold on number of channels';
sum_neg_threshold.tag     = 'sum_neg_threshold';
sum_neg_threshold.strtype = 'r';
sum_neg_threshold.num     = [1 1];
sum_neg_threshold.def     = @(val)nirs_get_defaults('preprocessNIRS.mark_negative.sum_neg_threshold', val{:});
sum_neg_threshold.help    = {['Enter the value of the threshold as ',...
    'a percentage of the total number of channels, so that ',...
    'data points which are marked bad for at least that many channels ',...
    'will be marked bad for all channels.']};

% Executable Branch
mark_negative      = cfg_exbranch;
mark_negative.name = 'Mark negative and interpolate';
mark_negative.tag  = 'mark_negative';
mark_negative.val  = {NIRSmat redo1 NIRSmatCopyChoice sum_neg_threshold};
mark_negative.prog = @nirs_run_mark_negative;
mark_negative.vout = @nirs_cfg_vout_mark_negative;
mark_negative.help = {['Preprocessing step: mark negative data points as ',...
    'bad data points, and interpolate from nearby values.']};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_mark_negative(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});