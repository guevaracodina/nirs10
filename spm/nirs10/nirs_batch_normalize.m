function nirs_batch_normalize(anatT1,dir_coreg)
[dirT1, fil, ext] = fileparts(anatT1);
tmp_file = fullfile(dirT1,['m' fil ext]);
if exist(tmp_file,'file')
    src_file = tmp_file;
else
    src_file = anatT1;
end
template = fullfile(spm('dir'),'templates','T1.nii');
%fwT1 = fullfile(dirT1,['w' fil ext(1:4)]);
fc1 =  fullfile(dirT1,['c1' fil ext(1:4)]);
if spm_existfile(fc1)
    c1_present = 1;
else
    c1_present = 0;
end
resample_files = {anatT1};
if c1_present
    resample_files = [resample_files; fc1];
end
owd = pwd;
cd(dir_coreg);
%Various options that we don't make available to the user in the GUI
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.source = {src_file};
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.wtsrc = '';
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = resample_files;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.template = {template};
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.smosrc = 8;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.smoref = 0;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.regtype = 'mni';
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.cutoff = 25;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.nits = 16;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = 1;
matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.preserve = 0;
matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.bb = [-78 -112 -50
    78 76 85];
matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.vox = [2 2 2];
matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.interp = 1;
matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.prefix = 'w';
spm_jobman('run',matlabbatch);
cd(owd);