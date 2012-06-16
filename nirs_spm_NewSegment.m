function nirs_spm_NewSegment(f)    
% Run SPM new segment module to segment anatomical image in different 
% tissues

%Options:
%Save corrected image only
save_corrected = [0 1];
%Native space
native_space = [1 0]; %for Dartel: [1 1]
%Warped tissue: none
warped_tissue = [0 0]; %for VBM (and Dartel?): [1 1]

[DirSPM,dummy,dummy2] = fileparts(which('spm')); 

% Go to anatomical image directory
[dir1 fil ext] = fileparts(f);
cd(dir1);

% Check if was run already
if spm_existfile(['c1' fil ext(1:end)]) || ...
        spm_existfile(['c1' fil ext(1:end-2)])
    disp(['In folder ' dir1 ', New Segment already run -- skipping']);    

else
    
    try
        clear matlabbatch
        matlabbatch{1}.spm.tools.preproc8.channel.vols = {f};
        % Default values
        matlabbatch{1}.spm.tools.preproc8.channel.biasreg = 0.0001;
        matlabbatch{1}.spm.tools.preproc8.channel.biasfwhm = 60;
        % Option to save bias field and bias corrected image
        matlabbatch{1}.spm.tools.preproc8.channel.write = save_corrected;
        
        % Gray matter
        matlabbatch{1}.spm.tools.preproc8.tissue(1).tpm = {fullfile(DirSPM, 'toolbox\Seg\TPM.nii,1')};
        % Number of gaussians used to represent the intensity distribution
        % of the tissue
        matlabbatch{1}.spm.tools.preproc8.tissue(1).ngaus = 2;
        % Option to produce both a tissue class image in alignment with the
        % otiginal (c*) and a one that can be used with the DARTEL toolbox
        % (rc*)
        matlabbatch{1}.spm.tools.preproc8.tissue(1).native = native_space;
        % Option to produce spatially normalised versions of the tissue 
        % class - both with (mwc*) and without (wc*) modulation
        matlabbatch{1}.spm.tools.preproc8.tissue(1).warped = warped_tissue;
        
        % Whitter matter
        matlabbatch{1}.spm.tools.preproc8.tissue(2).tpm = {fullfile(DirSPM, 'toolbox','Seg','TPM.nii,2')};
        matlabbatch{1}.spm.tools.preproc8.tissue(2).ngaus = 2;
        matlabbatch{1}.spm.tools.preproc8.tissue(2).native = native_space;
        matlabbatch{1}.spm.tools.preproc8.tissue(2).warped = warped_tissue;
        % CSF
        matlabbatch{1}.spm.tools.preproc8.tissue(3).tpm = {fullfile(DirSPM, 'toolbox','Seg','TPM.nii,3')};
        matlabbatch{1}.spm.tools.preproc8.tissue(3).ngaus = 2;
        matlabbatch{1}.spm.tools.preproc8.tissue(3).native = native_space;
        matlabbatch{1}.spm.tools.preproc8.tissue(3).warped = warped_tissue;
        % Bone
        matlabbatch{1}.spm.tools.preproc8.tissue(4).tpm = {fullfile(DirSPM, 'toolbox','Seg','TPM.nii,4')};
        matlabbatch{1}.spm.tools.preproc8.tissue(4).ngaus = 3;
        matlabbatch{1}.spm.tools.preproc8.tissue(4).native = native_space;
        matlabbatch{1}.spm.tools.preproc8.tissue(4).warped = warped_tissue;
        % Soft tissue
        matlabbatch{1}.spm.tools.preproc8.tissue(5).tpm = {fullfile(DirSPM, 'toolbox','Seg','TPM.nii,5')};
        matlabbatch{1}.spm.tools.preproc8.tissue(5).ngaus = 4;
        matlabbatch{1}.spm.tools.preproc8.tissue(5).native = native_space;
        matlabbatch{1}.spm.tools.preproc8.tissue(5).warped = warped_tissue;
        % Air/background
        matlabbatch{1}.spm.tools.preproc8.tissue(6).tpm = {fullfile(DirSPM, 'toolbox','Seg','TPM.nii,6')};
        matlabbatch{1}.spm.tools.preproc8.tissue(6).ngaus = 2;
        matlabbatch{1}.spm.tools.preproc8.tissue(6).native = native_space;
        matlabbatch{1}.spm.tools.preproc8.tissue(6).warped = warped_tissue;
        
        % Default options
        matlabbatch{1}.spm.tools.preproc8.warp.reg = 4;
        % (European brains)
        matlabbatch{1}.spm.tools.preproc8.warp.affreg = 'mni';
        matlabbatch{1}.spm.tools.preproc8.warp.samp = 3;
        % Option to write both forward (MNI->individual) and inverse
        % deformation fields
        % (written as .nii files with 3 volumes for x, y and z)
        matlabbatch{1}.spm.tools.preproc8.warp.write = warped_tissue;
        
        spm_jobman('run',matlabbatch);
    
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
        
        disp(['In folder ' dir1 ', New Segment failed to run.']);
    
    end
end


