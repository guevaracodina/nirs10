function nirs_batch_coreg(NIRS,newNIRSlocation)
if isfield(NIRS,'jobCoreg')
    job = NIRS.jobCoreg;
    job.force_redo = 1;
    job.NIRSmatCopyChoice = rmfield(job.NIRSmatCopyChoice,'NIRSmatCopy');
    job.NIRSmatCopyChoice.NIRSmatOverwrite = struct([]);
    nirs_run_coreg_new(job);
else
    clear matlabbatch
    matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.NIRSmat = {newNIRSlocation};
    matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.force_redo = 1;
    matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.NIRSmatCopyChoice.NIRSmatOverwrite = struct([]);
    matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.template_mode = 0;
    matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.anatT1 = {''};
    matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.segT1_4fit = {''};
    matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.anatT1_template = {'W:\spm8\templates\T1.nii'};
    matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.fid_in_subject_MNI = 0;
    matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.nasion_wMNI = [0 84 -48];
    matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.AL_wMNI = [-83 -19 -38];
    matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.AR_wMNI = [83 -19 -38];
    matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.GenDataTopo = 1;
    matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.render_choice.render_template = struct([]);
    matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.View6Projections = 0;
    matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.Save6Projections = 1;
    matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.ForceReprocess = 0;
    spm_jobman('run',matlabbatch);
end