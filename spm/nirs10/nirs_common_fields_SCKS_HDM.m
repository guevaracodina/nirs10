function [IOImat ROImat redo1 IOImatCopyChoice IC ...
    PhysioModel_Choice ROI_choice session_choice baseline_choice ...
    lpf_choice hpf_butter] = nirs_common_fields_SCKS_HDM(dir_name)

IOImat = ioi_dfg_IOImat(1);
redo1 = ioi_dfg_redo(0);
ROImat = ioi_dfg_ROImat(1);
IOImatCopyChoice = ioi_dfg_IOImatCopyChoice(dir_name);
IC = ioi_dfg_include_colors(0,1,1,1,1);

PhysioModel_Choice        = cfg_menu;
PhysioModel_Choice.name   = 'Choose hemodynamic model';
PhysioModel_Choice.tag    = 'PhysioModel_Choice';
PhysioModel_Choice.labels = {'Buxton-Friston','Zheng-Mayhew','Boas-Huppert'};
PhysioModel_Choice.values = {0,1,2};
PhysioModel_Choice.val    = {0};
PhysioModel_Choice.help   = {'Choose hemodynamic model'}';

ROI_choice = ioi_dfg_ROI_choice;
session_choice = ioi_dfg_session_choice;

hpf_butter = ioi_dfg_hpf_butter(1,0.01,3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
no_baseline_correction         = cfg_branch;
no_baseline_correction.tag     = 'no_baseline_correction';
no_baseline_correction.name    = 'No baseline correction';
no_baseline_correction.val     = {};
no_baseline_correction.help    = {'No correction will be made to baseline'};

baseline_percentile_HbR      = cfg_entry;
baseline_percentile_HbR.tag  = 'baseline_percentile_HbR';
baseline_percentile_HbR.name = 'Choose data percentile to set baseline for HbR';
baseline_percentile_HbR.strtype  = 'r';
baseline_percentile_HbR.num = [Inf Inf];
baseline_percentile_HbR.val{1} = 90;
baseline_percentile_HbR.help = {'Enter percentile of data (after filtering) to set baseline at, for HbR.'
    'Enter either a single number,'
    'to be applied to each selected session and ROI, or '
    'enter an array of offsets, by selected session and by selected ROI.'}';

baseline_percentile_HbT      = cfg_entry;
baseline_percentile_HbT.tag  = 'baseline_percentile_HbT';
baseline_percentile_HbT.name = 'Choose data percentile to set baseline for HbT';
baseline_percentile_HbT.strtype  = 'r';
baseline_percentile_HbT.num = [Inf Inf];
baseline_percentile_HbT.val{1} = 10;
baseline_percentile_HbT.help = {'Enter percentile of data (after filtering) to set baseline at, for HbT.'
    'Enter either a single number,'
    'to be applied to each selected session and ROI, or '
    'enter an array of offsets, by selected session and by selected ROI.'}';

baseline_percentile_flow      = cfg_entry;
baseline_percentile_flow.tag  = 'baseline_percentile_flow';
baseline_percentile_flow.name = 'Choose data percentile to set baseline for flow';
baseline_percentile_flow.strtype  = 'r';
baseline_percentile_flow.num = [Inf Inf];
baseline_percentile_flow.val{1} = 10;
baseline_percentile_flow.help = {'Enter percentile of data (after filtering) to set baseline at, for flow.'
    'Enter either a single number,'
    'to be applied to each selected session and ROI, or '
    'enter an array of offsets, by selected session and by selected ROI.'}';

baseline_percentile_choice         = cfg_branch;
baseline_percentile_choice.tag     = 'baseline_percentile_choice';
baseline_percentile_choice.name    = 'Baseline percentile choice';
baseline_percentile_choice.val     = {baseline_percentile_HbR baseline_percentile_HbT baseline_percentile_flow};
baseline_percentile_choice.help    = {'Set baseline to a chosen percentile'}';

baseline_offset_HbR      = cfg_entry;
baseline_offset_HbR.tag  = 'baseline_offset_HbR';
baseline_offset_HbR.name = 'Baseline offset for HbR';
baseline_offset_HbR.strtype  = 'r';
baseline_offset_HbR.num = [Inf Inf];
baseline_offset_HbR.val{1} = 0;
baseline_offset_HbR.help = {'Enter baseline offset for HbR. Enter either a single number,'
    'to be applied to each selected session and ROI, or  '
    'enter an array of offsets, by selected session and by selected ROI.'}';

baseline_offset_HbT      = cfg_entry;
baseline_offset_HbT.tag  = 'baseline_offset_HbT';
baseline_offset_HbT.name = 'Baseline offset for HbT';
baseline_offset_HbT.strtype  = 'r';
baseline_offset_HbT.num = [Inf Inf];
baseline_offset_HbT.val{1} = 0;
baseline_offset_HbT.help = {'Enter baseline offset for HbT. Enter either a single number,'
    'to be applied to each selected session and ROI, or  '
    'enter an array of offsets, by selected session and by selected ROI.'}';

baseline_offset_flow      = cfg_entry;
baseline_offset_flow.tag  = 'baseline_offset_flow';
baseline_offset_flow.name = 'Baseline offset for Flow';
baseline_offset_flow.strtype  = 'r';
baseline_offset_flow.num = [Inf Inf];
baseline_offset_flow.val{1} = 0;
baseline_offset_flow.help = {'Enter baseline offset for flow. Enter either a single number,'
    'to be applied to each selected session and ROI, or  '
    'enter an array of offsets, by selected session and by selected ROI.'}';

baseline_offset_choice         = cfg_branch;
baseline_offset_choice.tag     = 'baseline_offset_choice';
baseline_offset_choice.name    = 'Baseline offset choice';
baseline_offset_choice.val     = {baseline_offset_HbR baseline_offset_HbT baseline_offset_flow};
baseline_offset_choice.help    = {'Offset baseline by a chosen amount'};

baseline_choice        = cfg_choice;
baseline_choice.name   = 'Choose baseline selection method';
baseline_choice.tag    = 'baseline_choice';
baseline_choice.values = {no_baseline_correction,baseline_percentile_choice,baseline_offset_choice};
baseline_choice.val    = {baseline_percentile_choice};
baseline_choice.help   = {'Choose baseline selection method'}';
