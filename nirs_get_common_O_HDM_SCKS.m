function O = nirs_get_common_O_HDM_SCKS(job)
[O.all_sessions O.selected_sessions] = nirs_get_sessions(job);
[O.all_channels O.selected_channels] = nirs_get_channels(job);

%filters
O.HPF = nirs_get_HPF(job);
O.LPF = nirs_get_LPF(job);
O.use_onset_amplitudes = job.O.use_onset_amplitudes;
%Hemodynamic model choice
O.PhysioModel_Choice = job.O.Model_Choice; %eventually, Model_Choice should be
%replaced by PhysioModel_Choice in all the code.
O.prior_choice = job.O.prior_choice;
%Baseline choice
if isfield(job.O.baseline_choice,'baseline_percentile_choice')
    O.baseline_choice = 1;
    O.baseline_correction{1} = job.O.baseline_choice.baseline_percentile_choice.baseline_percentile_HbR/100;
    O.baseline_correction{2} = job.O.baseline_choice.baseline_percentile_choice.baseline_percentile_HbT/100;
else
    if isfield(job.O.baseline_choice,'baseline_offset_choice')
        O.baseline_choice = 2;
        O.baseline_correction{1} = job.O.baseline_choice.baseline_offset_choice.baseline_offset_HbR;
        O.baseline_correction{2} = job.O.baseline_choice.baseline_offset_choice.baseline_offset_HbT;
    else
        O.baseline_choice = 0;
        O.baseline_correction = [];
    end
end
O.baseline_HbR = job.O.baseline_HbR;
O.baseline_HbO = job.O.baseline_HbO;