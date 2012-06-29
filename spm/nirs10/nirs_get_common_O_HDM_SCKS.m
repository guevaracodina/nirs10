function O = nirs_get_common_O_HDM_SCKS(job)
[O.all_sessions O.selected_sessions] = nirs_get_sessions(job);
[O.all_channels O.selected_channels] = nirs_get_channels(job);

%filters
O.HPF = nirs_get_HPF(job);
O.LPF = nirs_get_LPF(job);

%Hemodynamic model choice
O.PhysioModel_Choice = job.PhysioModel_Choice;
%modalities to include
O.includeHbR = job.includeHbR;
O.includeHbO = job.includeHbO;
O.includeHbT = job.includeHbT;
%Baseline choice
if isfield(job.baseline_choice,'baseline_percentile_choice')
    O.baseline_choice = 1;
    O.baseline_correction{1} = job.baseline_choice.baseline_percentile_choice.baseline_percentile_HbR/100;
    O.baseline_correction{2} = job.baseline_choice.baseline_percentile_choice.baseline_percentile_HbT/100;
    O.baseline_correction{3} = job.baseline_choice.baseline_percentile_choice.baseline_percentile_flow/100;
else
    if isfield(job.baseline_choice,'baseline_offset_choice')
        O.baseline_choice = 2;
        O.baseline_correction{1} = job.baseline_choice.baseline_offset_choice.baseline_offset_HbR;
        O.baseline_correction{2} = job.baseline_choice.baseline_offset_choice.baseline_offset_HbT;
        O.baseline_correction{3} = job.baseline_choice.baseline_offset_choice.baseline_offset_flow;
    else
        O.baseline_choice = 0;
        O.baseline_correction = [];
    end
end