function HDM = nirs_run_HDM_cfg
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'HDM');
PhysioModel_Choice = nirs_dfg_Model_Choice;

[ IC
    PhysioModel_Choice ROI_choice session_choice baseline_choice ...
    lpf_choice hpf_butter] = ioi_common_fields_SCKS_HDM('HDM');

session_choice = nirs_dfg_session_choice;
channel_choice = nirs_dfg_channel_choice;
which_condition = nirs_dfg_which_condition;
%EM parameters
EM_parameters = nirs_dfg_hdm_EM;
[generate_figures save_figures] = nirs_dfg_generate_figures;
%Simulation options
simuOn = nirs_dfg_hdm_simu_options(0);
hdm_display_options = nirs_dfg_hdm_display_options;

use_onset_amplitudes      = cfg_menu;
use_onset_amplitudes.tag  = 'use_onset_amplitudes';
use_onset_amplitudes.name = 'Use onset amplitudes to weigh the hemodynamic response';
use_onset_amplitudes.labels = {'Yes','No'};
use_onset_amplitudes.values = {1,0};
use_onset_amplitudes.val  = {0};
use_onset_amplitudes.help = {'Use onset amplitudes as parameters to weigh the hemodynamic response.'}';

% Executable Branch
HDM      = cfg_exbranch;
HDM.name = 'Hemodynamic Modelling';
HDM.tag  = 'HDM';
HDM.val  = {NIRSmat redo1 NIRSmatCopyChoice PhysioModel_Choice ...
    which_condition session_choice channel_choice use_onset_amplitudes ...
    generate_figures save_figures hdm_display_options EM_parameters simuOn};
HDM.prog = @nirs_run_HDM;
HDM.vout = @nirs_cfg_vout_HDM;
HDM.help = {'Hemodynamic Modeling.'};

function vout = nirs_cfg_vout_HDM(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});