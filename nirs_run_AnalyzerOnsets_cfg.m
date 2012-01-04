function AnalyzerOnsets = nirs_run_AnalyzerOnsets_cfg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration  Read NIRS onsets to generate input to General Linear Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration  Read NIRS onsets for epilepsy from Analyzer 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

raw_onset_files        = cfg_files;
raw_onset_files.name    = 'Select onset files'; % The displayed name
raw_onset_files.tag     = 'raw_onset_files';
raw_onset_files.num     = [0 Inf];     % Number of inputs required
raw_onset_files.val{1}  = {''};
raw_onset_files.help    = {'Optional: Select raw onset files. '
    'Can be added at a later stage.'
    'Must specify one file for each data file, in same order.'}'; % help text displayed

freq_NIRS1      = cfg_entry;
freq_NIRS1.tag  = 'freq_NIRS1';
freq_NIRS1.name = 'Frequency of NIRS data for GLM';
freq_NIRS1.val{1} = []; %{1.9531}; %{19.5312};
freq_NIRS1.strtype = 'r';
freq_NIRS1.num     = [0 Inf];
freq_NIRS1.help    = {'Specify frequency of NIRS data (optional).'};

dp_NIRS1      = cfg_entry;
dp_NIRS1.tag  = 'dp_NIRS1';
dp_NIRS1.name = 'Number of data points in NIRS data for GLM';
dp_NIRS1.val{1} = []; %{1758};
dp_NIRS1.strtype = 'r';
dp_NIRS1.num     = [0 Inf];
dp_NIRS1.help    = {'Specify number of data time points in NIRS data for the GLM (optional).'};

NIRSmat_optional         = cfg_files; %Select NIRS.mat for this subject
NIRSmat_optional.name    = 'NIRS.mat'; % The displayed name
NIRSmat_optional.tag     = 'NIRSmat_optional';       %file names
NIRSmat_optional.filter  = 'mat';
NIRSmat_optional.ufilter = '^NIRS.mat$';
NIRSmat_optional.val{1}  = {''};
NIRSmat_optional.num     = [0 Inf];     % Number of inputs required
NIRSmat_optional.help    = {'Select NIRS.mat for the subject(s).'
    'If selecting more than one NIRS.mat, the onsets must have been already '
    'specified at an earlier stage. Otherwise, only one subject can be processed.'}'; % help text displayed


% Executable Branch
AnalyzerOnsets      = cfg_exbranch;
AnalyzerOnsets.name = 'Read NIRS onsets';
AnalyzerOnsets.tag  = 'AnalyzerOnsets';
AnalyzerOnsets.val  = {NIRSmat_optional raw_onset_files freq_NIRS1 dp_NIRS1};
AnalyzerOnsets.prog = @nirs_run_AnalyzerOnsets;
AnalyzerOnsets.vout = @nirs_cfg_vout_AnalyzerOnsets;
AnalyzerOnsets.help = {'Select NIRS structures (optional) and/or '
    'Analyzer 2 export files of onsets (also optional), '
    'to generate files of onsets and of confound regressors (pulse). '
    'This module can now be run by itself or as part of a larger batch.'
    'If specifying onset files, only one subject should be run.'}';

function vout = nirs_cfg_vout_AnalyzerOnsets(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});