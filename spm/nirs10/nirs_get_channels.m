function [all_channels selected_channels] = nirs_get_channels(job)
if isfield(job.channel_choice,'select_channels')
    all_channels = 0;
    selected_channels = job.channel_choice.select_channels.selected_channels;
else
    all_channels = 1;
    selected_channels = [];
end
