function AnalyzerOnsets = nirs_run_AnalyzerOnsets_cfg
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(0,'ons');

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

%**************************************************************************************************
%Ke Peng, include the heart rate repair option
%17/04/2012

avg_number         = cfg_entry;
avg_number.name    = 'Number of pulses to generate mean: ';
avg_number.tag     = 'avg_number';
avg_number.strtype = 'r';
avg_number.val     = {5};
avg_number.num     = [1 1];
%avg_number.def     = @(val)nirs_get_defaults(...
    %'model_specify.wls_bglm_specify.hpf_butter.hpf_butter_On.hpf_butter_freq', val{:});
avg_number.help    = {'Enter the number of pulses to generate the mean value for verification'
                      'This value has to be an integer.'}';

gap_def         = cfg_entry;
gap_def.name    = 'Times of mean value to define a gap: ';
gap_def.tag     = 'gap_def';
gap_def.strtype = 'r';
gap_def.val     = {1.8};
gap_def.num     = [1 1];
%avg_number.def     = @(val)nirs_get_defaults(...
    %'model_specify.wls_bglm_specify.hpf_butter.hpf_butter_On.hpf_butter_freq', val{:});
gap_def.help    = {'Enter the times of mean value to define a single gap.'};

cardiac_repair_on         = cfg_branch;
cardiac_repair_on.tag     = 'cardiac_repair_on';
cardiac_repair_on.name    = 'cardiac_repair_on';
cardiac_repair_on.val     = {avg_number gap_def};
cardiac_repair_on.help    = {'Please specify the parameters for gap filling.'};

cardiac_repair_off         = cfg_branch;
cardiac_repair_off.tag     = 'cardiac_repair_off';
cardiac_repair_off.name    = 'Cardiac onset repair off';
cardiac_repair_off.val     = {};
cardiac_repair_off.help    = {'Cardiac onset gap filling has been turned off.'};

cardiac_repair      = cfg_choice;
cardiac_repair.tag  = 'cardiac_repair';
cardiac_repair.name = 'Gap filling for Cadiac onsets (optional).';
%cardiac_repair.labels = {'Yes','No'};
cardiac_repair.values = {cardiac_repair_on cardiac_repair_off};
cardiac_repair.val = {cardiac_repair_on};
cardiac_repair.help = {'Choose whether to fill the gaps of the cardiac onsets is needed.'
                       'This process compares each value in cardiac onsets with the mean of certain previous cardiac interval values'
                       'interpolate values to fill the gaps between two cardiac pulses'}';

onset_to_keep         = cfg_entry;
onset_to_keep.name    = 'Label for single onset type to keep';
onset_to_keep.tag     = 'onset_to_keep';
onset_to_keep.strtype = 's';
onset_to_keep.val     = {''};
onset_to_keep.num     = [1 Inf];
onset_to_keep.help    = {'Enter empty string to process as usual (i.e. keep onsets)'
    'Alternatively, enter the name of the onset type (e..g spkLT ) in order'
    'to exclude other onsets and keep only this onset type.'}';

%*****************************************************************************************************
% Executable Branch
AnalyzerOnsets      = cfg_exbranch;
AnalyzerOnsets.name = 'Read NIRS onsets';
AnalyzerOnsets.tag  = 'AnalyzerOnsets';
AnalyzerOnsets.val  = {NIRSmat redo1 NIRSmatCopyChoice raw_onset_files ...
    onset_to_keep freq_NIRS1 dp_NIRS1 cardiac_repair};
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