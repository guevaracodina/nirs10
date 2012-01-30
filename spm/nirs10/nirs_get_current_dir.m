function [newDir job] = nirs_get_current_dir(job,SubjIdx)
%returns current directory and updated location of NIRS.mat
[dir_NIRSmat dummy] = fileparts(job.NIRSmat{SubjIdx});
if isfield(job.NIRSmatCopyChoice,'NIRSmatCopy')
    newDir = job.NIRSmatCopyChoice.NIRSmatCopy.NewNIRSdir;
    newDir = fullfile(dir_NIRSmat,newDir);
    if ~exist(newDir,'dir'),mkdir(newDir); end
    job.NIRSmat{SubjIdx} = fullfile(newDir,'NIRS.mat');
else
    newDir = dir_NIRSmat;
end