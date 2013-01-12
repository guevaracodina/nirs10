function target_sampling_rate = nirs_dfg_target_sampling_rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sampling rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampling_rate      = cfg_entry;
sampling_rate.tag  = 'sampling_rate';
sampling_rate.name = 'Sampling rate in Hz for the purpose of TrRVRV and erdf calculation';
sampling_rate.val  = {1};
sampling_rate.strtype = 'r';
sampling_rate.num     = [1 1];
sampling_rate.help    = {'Specify sampling rate in Hz.'
    'This is used only for the purpose of TrRVRV and erdf calculation'
    'This option is used to speed up estimation of the GLM.'}';

specified_sampling_rate         = cfg_branch;
specified_sampling_rate.tag     = 'specified_sampling_rate';
specified_sampling_rate.name    = 'Specify sampling rate';
specified_sampling_rate.val     = {sampling_rate};
specified_sampling_rate.help    = {'Specified sampling rate'};

raw_data_sampling_rate         = cfg_branch;
raw_data_sampling_rate.tag     = 'raw_data_sampling_rate';
raw_data_sampling_rate.name    = 'Raw data sampling rate';
raw_data_sampling_rate.val     = {};
raw_data_sampling_rate.help    = {};

target_sampling_rate        = cfg_choice;
target_sampling_rate.name   = 'Choose sampling rate for the purpose of TrRVRV and erdf calculation';
target_sampling_rate.tag    = 'target_sampling_rate';
target_sampling_rate.values = {raw_data_sampling_rate,specified_sampling_rate};
target_sampling_rate.val    = {raw_data_sampling_rate};
target_sampling_rate.help   = {'Choose whether to specify a sampling rate for the data in the GLM.'
    'for the purpose of TrRVRV and erdf calculation.'
    'It will not affect the actual sampling rate'}';
