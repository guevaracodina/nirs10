function HDM = nirs_run_HDM_cfg
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'HDM');
O = nirs_common_fields_SCKS_HDM;
IC = nirs_dfg_include_colors(1,1,1); %colors to include
session_choice = nirs_dfg_session_choice;
channel_choice = nirs_dfg_channel_choice;
which_condition = nirs_dfg_which_condition;
EM_parameters = nirs_dfg_hdm_EM; %EM parameters
simuOn = nirs_dfg_hdm_simu_options(0); %Simulation options
hdm_display_options = nirs_dfg_hdm_display_options;
lpf_choice = nirs_dfg_lpf_choice(1,1.5);
hpf_butter = nirs_dfg_hpf_butter(1,0.01,3);

% Executable Branch
HDM      = cfg_exbranch;
HDM.name = 'Hemodynamic Modelling';
HDM.tag  = 'HDM';
HDM.val  = {NIRSmat redo1 NIRSmatCopyChoice O ...
    which_condition session_choice channel_choice lpf_choice hpf_butter ...    
    IC hdm_display_options EM_parameters simuOn};
HDM.prog = @nirs_run_HDM;
HDM.vout = @nirs_cfg_vout_HDM;
HDM.help = {'Hemodynamic Modeling.'};

function vout = nirs_cfg_vout_HDM(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});