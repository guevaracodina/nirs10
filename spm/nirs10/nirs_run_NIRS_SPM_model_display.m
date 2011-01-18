function out = nirs_run_NIRS_SPM_model_display(job)
%Load NIRS.mat information
clear NIRS
load(job.NIRSmat{1,1});


save(job.NIRSmat{1,1},'NIRS');
out.NIRSmat{1} = fullfile(NIRS.subj_path,'NIRS.mat');
end