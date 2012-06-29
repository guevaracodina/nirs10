function [all_sessions selected_sessions] = nirs_get_sessions(job)
if isfield(job.session_choice,'select_sessions')
    all_sessions = 0;
    selected_sessions = job.session_choice.select_sessions.selected_sessions;
else
    all_sessions = 1;
    selected_sessions = [];
end
