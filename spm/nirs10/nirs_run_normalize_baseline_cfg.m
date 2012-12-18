function normalize_baseline = nirs_run_normalize_baseline_cfg
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'base');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.4. Normalize to baseline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Choose to normalize Optical Densities to initial value or to median
Normalize_OD = cfg_menu;
Normalize_OD.tag  = 'Normalize_OD';
Normalize_OD.name = 'Normalization method';
Normalize_OD.labels = {'Median','Last Value','Mean',...
    'Median- use on session 1 for all','Mean of first 60 seconds only',...
    'Mean of first 10 seconds only','Mean of first 60 seconds only from the first session only (used for all other sessions)'};
Normalize_OD.values = {0,1,2,3,4,5,6};
Normalize_OD.def  = @(val)nirs_get_defaults('preprocessNIRS.normalize_baseline.Normalize_OD', val{:});
Normalize_OD.help = {'Choose normalization of Optical Densities'
    'Median (preferred), Initial value, or Mean.'
    'Median - all sessions: this will normalize each channel separately by the'
    'value of the median of the first session'}';

add_or_mult      = cfg_menu;
add_or_mult.tag  = 'add_or_mult';
add_or_mult.name = 'Additive or multiplicative';
add_or_mult.labels = {'Additive', 'Multiplicative'};
add_or_mult.values = {1,0};
add_or_mult.def  = @(val)nirs_get_defaults('preprocessNIRS.normalize_baseline.add_or_mult', val{:});
add_or_mult.help = {'Select whether using additive (on concentrations for example)'
    'or multiplicative (on optical intensities) normalization.' }';

baseline_duration         = cfg_entry;
baseline_duration.name    = 'Baseline duration'; % The displayed name
baseline_duration.tag     = 'baseline_duration';       %file names
baseline_duration.strtype = 'r';
baseline_duration.num     = [1 1];     % Number of inputs required
%baseline_duration.val{1}  = 100;
baseline_duration.def     = @(val)nirs_get_defaults('preprocessNIRS.normalize_baseline.baseline_duration', val{:});
baseline_duration.help    = {'Enter the baseline duration in seconds to use '
    'prior to stimuli - applies only '
    'to the normalization type by stimuli below.)'}';

normalization_type      = cfg_menu;
normalization_type.tag  = 'normalization_type';
normalization_type.name = 'Normalization type';
normalization_type.labels = {'Global', 'By bad point segments', 'By stimuli'};
normalization_type.values = {1,2,3};
normalization_type.def  = @(val)nirs_get_defaults('preprocessNIRS.normalize_baseline.normalization_type', val{:});
normalization_type.help = {'Normalization type: global, by bad point segments, before stimuli.'
    'When normalizing by stimuli, the code finds each stimulus instance '
    'and uses the time window specified prior to the stimulus for normalization.'}';

Analyzer_sf         = cfg_entry;
Analyzer_sf.name    = 'Scaling factor'; % The displayed name
Analyzer_sf.tag     = 'Analyzer_sf';       %file names
Analyzer_sf.strtype = 'r';
Analyzer_sf.num     = [1 1];     % Number of inputs required
%Analyzer_sf.val{1}  = 100;
Analyzer_sf.def     = @(val)nirs_get_defaults('preprocessNIRS.normalize_baseline.Analyzer_sf', val{:});
Analyzer_sf.help    = {'Apply a scaling factor on the amplitude of '
    'all channels (for easier visualization with Analyzer.)'}';

nirs_filling_jumps = nirs_dfg_nirs_filling_jumps;

% Executable Branch
normalize_baseline      = cfg_exbranch;
normalize_baseline.name = 'Normalize Baseline';
normalize_baseline.tag  = 'normalize_baseline';
normalize_baseline.val  = {NIRSmat redo1 NIRSmatCopyChoice Normalize_OD add_or_mult ...
    baseline_duration normalization_type Analyzer_sf nirs_filling_jumps};
normalize_baseline.prog = @nirs_run_normalize_baseline;
normalize_baseline.vout = @nirs_cfg_vout_normalize_baseline;
normalize_baseline.help = {'Normalize to baseline'}';

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_normalize_baseline(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
