function O = nirs_common_fields_SCKS_HDM

PhysioModel_Choice = nirs_dfg_Model_Choice;

use_onset_amplitudes      = cfg_menu;
use_onset_amplitudes.tag  = 'use_onset_amplitudes';
use_onset_amplitudes.name = 'Use onset amplitudes to weigh the hemodynamic response';
use_onset_amplitudes.labels = {'Yes','No'};
use_onset_amplitudes.values = {1,0};
use_onset_amplitudes.val  = {0};
use_onset_amplitudes.help = {'Use onset amplitudes as parameters to weigh the hemodynamic response.'}';

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

baseline_percentile_choice         = cfg_branch;
baseline_percentile_choice.tag     = 'baseline_percentile_choice';
baseline_percentile_choice.name    = 'Baseline percentile choice';
baseline_percentile_choice.val     = {baseline_percentile_HbR baseline_percentile_HbT};
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

baseline_HbR      = cfg_entry;
baseline_HbR.tag  = 'baseline_HbR';
baseline_HbR.name = 'Baseline HbR concentration in micro molar';
baseline_HbR.strtype  = 'r';
baseline_HbR.num = [1 1];
baseline_HbR.val{1} = 40;
baseline_HbR.help = {''}';

baseline_HbO      = cfg_entry;
baseline_HbO.tag  = 'baseline_HbO';
baseline_HbO.name = 'Baseline HbO concentration in micro molar';
baseline_HbO.strtype  = 'r';
baseline_HbO.num = [1 1];
baseline_HbO.val{1} = 60;
baseline_HbO.help = {''}';

baseline_offset_choice         = cfg_branch;
baseline_offset_choice.tag     = 'baseline_offset_choice';
baseline_offset_choice.name    = 'Baseline offset choice';
baseline_offset_choice.val     = {baseline_offset_HbR baseline_offset_HbT};
baseline_offset_choice.help    = {'Offset baseline by a chosen amount'};

baseline_choice        = cfg_choice;
baseline_choice.name   = 'Choose baseline selection method';
baseline_choice.tag    = 'baseline_choice';
baseline_choice.values = {no_baseline_correction,baseline_percentile_choice,baseline_offset_choice};
baseline_choice.val    = {baseline_percentile_choice};
baseline_choice.help   = {'Choose baseline selection method'}';

prior_choice = nirs_dfg_prior_choice;

O         = cfg_branch;
O.name     = 'Options';
O.tag    = 'O';
O.val     = {PhysioModel_Choice prior_choice baseline_choice ...
    use_onset_amplitudes baseline_HbR baseline_HbO}; 
O.help    = {'Choose various options.'};