function nirs_seizure_evolution(OP)
pathbatch = OP.pathbatch;
try
    matlabbatch{1}.spm.tools.nirs10.model_specify.liom_intrasubject_average.NIRSmat{1} = OP.NIRSmat;
    matlabbatch{1}.spm.tools.nirs10.model_specify.liom_intrasubject_average.force_redo = OP.force_redo;
    matlabbatch{1}.spm.tools.nirs10.model_specify.liom_intrasubject_average.NIRSmatCopyChoice.NIRSmatCopy.NewNIRSdir = OP.sz_dir;
    matlabbatch{1}.spm.tools.nirs10.model_specify.liom_intrasubject_average.sessions = OP.session;
    matlabbatch{1}.spm.tools.nirs10.model_specify.liom_intrasubject_average.subj.input_onsets{1} = OP.onsets;
    matlabbatch{1}.spm.tools.nirs10.model_specify.liom_intrasubject_average.subj.multi_reg = {''};
    matlabbatch{1}.spm.tools.nirs10.model_specify.liom_intrasubject_average.units = 1;
    matlabbatch{1}.spm.tools.nirs10.model_specify.liom_intrasubject_average.time_res = 1;
    matlabbatch{1}.spm.tools.nirs10.model_specify.liom_intrasubject_average.GLM_include_cardiac = 0;
    matlabbatch{1}.spm.tools.nirs10.model_specify.liom_intrasubject_average.GLM_include_Mayer = 0;
    matlabbatch{1}.spm.tools.nirs10.model_specify.liom_intrasubject_average.vasomotion_choice.no_vasomotion = struct([]);
    matlabbatch{1}.spm.tools.nirs10.model_specify.liom_intrasubject_average.NIRSchannelsConfound.NoNIRSconfounds = struct([]);
    matlabbatch{1}.spm.tools.nirs10.model_specify.liom_intrasubject_average.GenerateHbT = 1;
    matlabbatch{1}.spm.tools.nirs10.model_specify.liom_intrasubject_average.channel_pca = 0;
    matlabbatch{1}.spm.tools.nirs10.model_specify.liom_intrasubject_average.hpf_butter.hpf_butter_On.hpf_butter_freq = 0.01;
    matlabbatch{1}.spm.tools.nirs10.model_specify.liom_intrasubject_average.hpf_butter.hpf_butter_On.hpf_butter_order = 2;
    matlabbatch{1}.spm.tools.nirs10.model_specify.liom_intrasubject_average.AvgFilters.nirs_hpf.hpf_none = struct([]);
    matlabbatch{1}.spm.tools.nirs10.model_specify.liom_intrasubject_average.AvgFilters.nirs_lpf.lpf_gauss.fwhm1 = 1.5;
    matlabbatch{1}.spm.tools.nirs10.model_specify.liom_intrasubject_average.averaging_choice.seizure_evolution.onset_delays = OP.delays;
    matlabbatch{1}.spm.tools.nirs10.model_specify.liom_intrasubject_average.averaging_choice.seizure_evolution.onset_duration = OP.onset_duration;
    switch OP.base_choice
        case 1
            matlabbatch{1}.spm.tools.nirs10.model_specify.liom_intrasubject_average.baseline_choice.baseline_block_averaging.baseline_offset = OP.base_offset;
            matlabbatch{1}.spm.tools.nirs10.model_specify.liom_intrasubject_average.baseline_choice.baseline_block_averaging.baseline_duration = OP.base_duration;
        case 3           
            matlabbatch{1}.spm.tools.nirs10.model_specify.liom_intrasubject_average.baseline_choice.unique_baseline.baseline_start = OP.baseline_start;
            matlabbatch{1}.spm.tools.nirs10.model_specify.liom_intrasubject_average.baseline_choice.unique_baseline.baseline_duration = OP.baseline_duration;            
            matlabbatch{1}.spm.tools.nirs10.model_specify.liom_intrasubject_average.baseline_choice.unique_baseline.baseline_session = 1; %OP.session;
    end
    matlabbatch{2}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1) = cfg_dep;
    matlabbatch{2}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).tname = 'NIRS.mat';
    matlabbatch{2}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).tgt_spec{1}(1).name = 'filter';
    matlabbatch{2}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).tgt_spec{1}(1).value = 'mat';
    matlabbatch{2}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{2}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{2}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).sname = 'LIOM Intrasubject Average: NIRS.mat';
    matlabbatch{2}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{2}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).src_output = substruct('.','NIRSmat');
    matlabbatch{2}.spm.tools.nirs10.model_estimate.liom_contrast.force_redo = 1;
    matlabbatch{2}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmatCopyChoice.NIRSmatOverwrite = struct([]);
    matlabbatch{2}.spm.tools.nirs10.model_estimate.liom_contrast.ContrastChoice.automated_contrasts.NonlinearEpilepsyOn = 0;
    matlabbatch{2}.spm.tools.nirs10.model_estimate.liom_contrast.Sessions = ''; %use the only session left after intrasubject_average
    matlabbatch{2}.spm.tools.nirs10.model_estimate.liom_contrast.view = OP.views;
    matlabbatch{2}.spm.tools.nirs10.model_estimate.liom_contrast.TopoData = {''};
    matlabbatch{2}.spm.tools.nirs10.model_estimate.liom_contrast.StatMethod = 0;
    matlabbatch{2}.spm.tools.nirs10.model_estimate.liom_contrast.UseCorrelRes = 0;
    matlabbatch{2}.spm.tools.nirs10.model_estimate.liom_contrast.GenerateStats = 1;
    matlabbatch{2}.spm.tools.nirs10.model_estimate.liom_contrast.AllowExtrapolation = 1;
    matlabbatch{2}.spm.tools.nirs10.model_estimate.liom_contrast.spatial_LPF.spatial_LPF_Off = struct([]);
    matlabbatch{2}.spm.tools.nirs10.model_estimate.liom_contrast.contrast_p_value = OP.p_value;
    matlabbatch{2}.spm.tools.nirs10.model_estimate.liom_contrast.display_options.GenerateInverted = 1;
    matlabbatch{2}.spm.tools.nirs10.model_estimate.liom_contrast.display_options.GroupColorbars = 0;
    matlabbatch{2}.spm.tools.nirs10.model_estimate.liom_contrast.display_options.contrast_figures = 3;
    matlabbatch{2}.spm.tools.nirs10.model_estimate.liom_contrast.display_options.override_colorbar.colorbar_default = struct([]);
    matlabbatch{2}.spm.tools.nirs10.model_estimate.liom_contrast.display_options.figures_visible = 0;
    matlabbatch{2}.spm.tools.nirs10.model_estimate.liom_contrast.display_options.GroupFiguresIntoSubplots = 1;
    matlabbatch{2}.spm.tools.nirs10.model_estimate.liom_contrast.display_options.output_unc = 1;
    matlabbatch{2}.spm.tools.nirs10.model_estimate.liom_contrast.display_options.SmallFigures = 1;
    matlabbatch{2}.spm.tools.nirs10.model_estimate.liom_contrast.display_options.write_neg_pos = 0;
    matlabbatch{2}.spm.tools.nirs10.model_estimate.liom_contrast.display_options.save_nifti_contrasts = 0;save(fullfile(pathbatch,'Last_batch_sz_evolution'),'matlabbatch');
    spm_jobman('run',matlabbatch);
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
end