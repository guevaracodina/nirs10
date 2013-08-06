function SCKS = nirs_run_SCKS_cfg
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'SCKS');
O = nirs_common_fields_SCKS_HDM;
IC = nirs_dfg_include_colors(1,1,1); %colors to include
session_choice = nirs_dfg_session_choice;
channel_choice = nirs_dfg_channel_choice;
which_condition = nirs_dfg_which_condition;
SCKS_parameters = nirs_dfg_SCKSparams; %SCKS parameters
simuOn = nirs_dfg_hdm_simu_options(0); %Simulation options
SCKS_display_options = nirs_dfg_SCKS_display_options;
lpf_choice = nirs_dfg_lpf_choice(0,1.5);
hpf_filter = nirs_dfg_hpf_filter;
target_sampling_rate = nirs_dfg_target_sampling_rate_HDM;

% Executable Branch
SCKS      = cfg_exbranch;
SCKS.name = 'SCKS estimation';
SCKS.tag  = 'SCKS';
SCKS.val  = {NIRSmat redo1 NIRSmatCopyChoice O ...
    which_condition session_choice channel_choice target_sampling_rate lpf_choice hpf_filter ...    
    IC SCKS_display_options SCKS_parameters simuOn};
SCKS.prog = @nirs_run_SCKS;
SCKS.vout = @nirs_cfg_vout_SCKS;
SCKS.help = {'SCKS Estimation.'};

function vout = nirs_cfg_vout_SCKS(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});