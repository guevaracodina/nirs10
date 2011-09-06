%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4252 $)
%-----------------------------------------------------------------------
matlabbatch{1}.spm.tools.nirs10.readNIRS.criugm1.study_cfg.study_path = '<UNDEFINED>';
matlabbatch{1}.spm.tools.nirs10.readNIRS.criugm1.study_cfg.indvdata.indvdata_chosen = struct([]);
matlabbatch{1}.spm.tools.nirs10.readNIRS.criugm1.subj(1).subj_id = '001';
matlabbatch{1}.spm.tools.nirs10.readNIRS.criugm1.subj(1).age1 = 25;
matlabbatch{1}.spm.tools.nirs10.readNIRS.criugm1.subj(1).anatT1 = {'W:\Claudine\RawData\Anatomiques_CIHR\CIHR_MDEIE_001\s201007051100-0002-00001-000160-01.nii'};
matlabbatch{1}.spm.tools.nirs10.readNIRS.criugm1.subj(1).helmet.text_brainsight = {'W:\Claudine\RawData\S001\S001_2010-07-07-09h30_optodesPos.txt  '};
matlabbatch{1}.spm.tools.nirs10.readNIRS.criugm1.subj(1).CWsystem = 6;
matlabbatch{1}.spm.tools.nirs10.readNIRS.criugm1.subj(1).nirs_files = {'W:\Claudine\RawData\S001\S001_StroopBlocs_1.nirs'};
matlabbatch{1}.spm.tools.nirs10.readNIRS.criugm1.subj(1).protocol = {'W:\Claudine\Protocol_ex_baseline.mat'};
matlabbatch{1}.spm.tools.nirs10.readNIRS.criugm1.subj(1).TopoData = {''};
matlabbatch{1}.spm.tools.nirs10.readNIRS.criugm1.subj(1).boldmask = {''};
matlabbatch{1}.spm.tools.nirs10.readNIRS.criugm1.subj(2).subj_id = '002';
matlabbatch{1}.spm.tools.nirs10.readNIRS.criugm1.subj(2).age1 = 25;
matlabbatch{1}.spm.tools.nirs10.readNIRS.criugm1.subj(2).anatT1 = {'W:\Claudine\RawData\Anatomiques_CIHR\CIHR_MDEIE_002\s201007090900-0002-00001-000160-01.nii'};
matlabbatch{1}.spm.tools.nirs10.readNIRS.criugm1.subj(2).helmet.text_brainsight = {'W:\Claudine\RawData\S002\S002_2010-07-09_12h00_optodesPos.txt  '};
matlabbatch{1}.spm.tools.nirs10.readNIRS.criugm1.subj(2).CWsystem = 6;
matlabbatch{1}.spm.tools.nirs10.readNIRS.criugm1.subj(2).nirs_files = {'W:\Claudine\RawData\S002\Mod_S002_StroopBlocs_1.nirs'};
matlabbatch{1}.spm.tools.nirs10.readNIRS.criugm1.subj(2).protocol = {'W:\Claudine\Protocol_ex_baseline.mat'};
matlabbatch{1}.spm.tools.nirs10.readNIRS.criugm1.subj(2).TopoData = {''};
matlabbatch{1}.spm.tools.nirs10.readNIRS.criugm1.subj(2).boldmask = {''};
n=2;
matlabbatch{n}.spm.tools.nirs10.preprocANAT.MCsegment1.NIRSmat_optional(1) = cfg_dep;
matlabbatch{n}.spm.tools.nirs10.preprocANAT.MCsegment1.NIRSmat_optional(1).tname = 'NIRS.mat';
matlabbatch{n}.spm.tools.nirs10.preprocANAT.MCsegment1.NIRSmat_optional(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{n}.spm.tools.nirs10.preprocANAT.MCsegment1.NIRSmat_optional(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{n}.spm.tools.nirs10.preprocANAT.MCsegment1.NIRSmat_optional(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{n}.spm.tools.nirs10.preprocANAT.MCsegment1.NIRSmat_optional(1).tgt_spec{1}(2).value = 'e';
matlabbatch{n}.spm.tools.nirs10.preprocANAT.MCsegment1.NIRSmat_optional(1).sname = 'Read and format CRIUGM data: NIRS.mat';
matlabbatch{n}.spm.tools.nirs10.preprocANAT.MCsegment1.NIRSmat_optional(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{n}.spm.tools.nirs10.preprocANAT.MCsegment1.NIRSmat_optional(1).src_output = substruct('.','NIRSmat');
matlabbatch{n}.spm.tools.nirs10.preprocANAT.MCsegment1.DelPreviousData = 0;
matlabbatch{n}.spm.tools.nirs10.preprocANAT.MCsegment1.NewDirCopyNIRS.CreateNIRSCopy_false = struct([]);
matlabbatch{n}.spm.tools.nirs10.preprocANAT.MCsegment1.image_in = {''}; %{Anatfile};
matlabbatch{n}.spm.tools.nirs10.preprocANAT.MCsegment1.output_autonaming = 0;
matlabbatch{n}.spm.tools.nirs10.preprocANAT.MCsegment1.output_prefix = 'Only edit if you chose ''No'' to ''Automatic output naming''';
matlabbatch{n}.spm.tools.nirs10.preprocANAT.MCsegment1.skn.sorting_method = 1;
matlabbatch{n}.spm.tools.nirs10.preprocANAT.MCsegment1.skl.sorting_method = 2;
matlabbatch{n}.spm.tools.nirs10.preprocANAT.MCsegment1.csf.sorting_method = 0;
matlabbatch{n}.spm.tools.nirs10.preprocANAT.MCsegment1.grm.sorting_method = 0;
matlabbatch{n}.spm.tools.nirs10.preprocANAT.MCsegment1.wtm.sorting_method = 0;
matlabbatch{n}.spm.tools.nirs10.preprocANAT.MCsegment1.thresh_as = 0.6;
matlabbatch{n}.spm.tools.nirs10.preprocANAT.MCsegment1.head_shadow.thresh_hs = 0.6;
matlabbatch{n}.spm.tools.nirs10.preprocANAT.MCsegment1.head_shadow.se_size_hs = 2;
matlabbatch{n}.spm.tools.nirs10.preprocANAT.MCsegment1.rebel_surrounding = 3;
matlabbatch{n}.spm.tools.nirs10.preprocANAT.MCsegment1.rebel_thresh_hs = 0.3;
matlabbatch{n}.spm.tools.nirs10.preprocANAT.MCsegment1.process_image.se_size_pi = 2;
matlabbatch{n}.spm.tools.nirs10.preprocANAT.MCsegment1.process_image.gaussfilt_size = 7;
matlabbatch{n}.spm.tools.nirs10.preprocANAT.MCsegment1.process_image.gaussfilt_sdev = 4;
%%%%%%%%%%
matlabbatch{3}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.NIRSmat(1) = cfg_dep;
matlabbatch{3}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.NIRSmat(1).tname = 'NIRS.mat';
matlabbatch{3}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.NIRSmat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{3}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.NIRSmat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{3}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.NIRSmat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{3}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.NIRSmat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{3}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.NIRSmat(1).sname = 'MC Segmentation: NIRS.mat';
matlabbatch{3}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.NIRSmat(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{3}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.NIRSmat(1).src_output = substruct('.','NIRSmat');
matlabbatch{3}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.DelPreviousData = 0;
matlabbatch{3}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.NewDirCopyNIRS.CreateNIRSCopy_false = struct([]);
matlabbatch{3}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.Normalize_OD = 0;
matlabbatch{3}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.add_or_mult = 0;
matlabbatch{3}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.baseline_duration = 2;
matlabbatch{3}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.normalization_type = 2;
matlabbatch{3}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.Analyzer_sf = 1;
matlabbatch{4}.spm.tools.nirs10.preprocessNIRS.ODtoHbOHbR.NIRSmat(1) = cfg_dep;
matlabbatch{4}.spm.tools.nirs10.preprocessNIRS.ODtoHbOHbR.NIRSmat(1).tname = 'NIRS.mat';
matlabbatch{4}.spm.tools.nirs10.preprocessNIRS.ODtoHbOHbR.NIRSmat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{4}.spm.tools.nirs10.preprocessNIRS.ODtoHbOHbR.NIRSmat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{4}.spm.tools.nirs10.preprocessNIRS.ODtoHbOHbR.NIRSmat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{4}.spm.tools.nirs10.preprocessNIRS.ODtoHbOHbR.NIRSmat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{4}.spm.tools.nirs10.preprocessNIRS.ODtoHbOHbR.NIRSmat(1).sname = 'Normalize Baseline: NIRS.mat';
matlabbatch{4}.spm.tools.nirs10.preprocessNIRS.ODtoHbOHbR.NIRSmat(1).src_exbranch = substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{4}.spm.tools.nirs10.preprocessNIRS.ODtoHbOHbR.NIRSmat(1).src_output = substruct('.','NIRSmat');
matlabbatch{4}.spm.tools.nirs10.preprocessNIRS.ODtoHbOHbR.DelPreviousData = 0;
matlabbatch{4}.spm.tools.nirs10.preprocessNIRS.ODtoHbOHbR.NewDirCopyNIRS.CreateNIRSCopy_false = struct([]);
matlabbatch{4}.spm.tools.nirs10.preprocessNIRS.ODtoHbOHbR.PVF = [50 50];
matlabbatch{4}.spm.tools.nirs10.preprocessNIRS.ODtoHbOHbR.nirs_lpf2.lpf_none = struct([]);
matlabbatch{5}.spm.tools.nirs10.coregNIRS.coreg1.NIRSmat(1) = cfg_dep;
matlabbatch{5}.spm.tools.nirs10.coregNIRS.coreg1.NIRSmat(1).tname = 'NIRS.mat';
matlabbatch{5}.spm.tools.nirs10.coregNIRS.coreg1.NIRSmat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{5}.spm.tools.nirs10.coregNIRS.coreg1.NIRSmat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{5}.spm.tools.nirs10.coregNIRS.coreg1.NIRSmat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{5}.spm.tools.nirs10.coregNIRS.coreg1.NIRSmat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{5}.spm.tools.nirs10.coregNIRS.coreg1.NIRSmat(1).sname = 'Convert OD to HbO/HbR : NIRS.mat';
matlabbatch{5}.spm.tools.nirs10.coregNIRS.coreg1.NIRSmat(1).src_exbranch = substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{5}.spm.tools.nirs10.coregNIRS.coreg1.NIRSmat(1).src_output = substruct('.','NIRSmat');
matlabbatch{5}.spm.tools.nirs10.coregNIRS.coreg1.DelPreviousData = 0;
matlabbatch{5}.spm.tools.nirs10.coregNIRS.coreg1.NewDirCopyNIRS.CreateNIRSCopy_false = struct([]);
matlabbatch{5}.spm.tools.nirs10.coregNIRS.coreg1.anatT1 = {''};
matlabbatch{5}.spm.tools.nirs10.coregNIRS.coreg1.segT1_4fit = {''};
matlabbatch{5}.spm.tools.nirs10.coregNIRS.coreg1.anatT1_template = {'D:\Users\Philippe Pouliot\spm8\templates\T1.nii'};
matlabbatch{5}.spm.tools.nirs10.coregNIRS.coreg1.fid_in_subject_MNI = 0;
matlabbatch{5}.spm.tools.nirs10.coregNIRS.coreg1.nasion_wMNI = [0 84 -48];
matlabbatch{5}.spm.tools.nirs10.coregNIRS.coreg1.AL_wMNI = [-83 -19 -38];
matlabbatch{5}.spm.tools.nirs10.coregNIRS.coreg1.AR_wMNI = [83 -19 -38];
matlabbatch{5}.spm.tools.nirs10.coregNIRS.coreg1.GenDataTopo = 1;
matlabbatch{5}.spm.tools.nirs10.coregNIRS.coreg1.render_choice.render_template = struct([]);
matlabbatch{5}.spm.tools.nirs10.coregNIRS.coreg1.View6Projections = 1;
matlabbatch{6}.spm.tools.nirs10.model_specify.wls_bglm_specify.NIRSmat(1) = cfg_dep;
matlabbatch{6}.spm.tools.nirs10.model_specify.wls_bglm_specify.NIRSmat(1).tname = 'NIRS.mat';
matlabbatch{6}.spm.tools.nirs10.model_specify.wls_bglm_specify.NIRSmat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{6}.spm.tools.nirs10.model_specify.wls_bglm_specify.NIRSmat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{6}.spm.tools.nirs10.model_specify.wls_bglm_specify.NIRSmat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{6}.spm.tools.nirs10.model_specify.wls_bglm_specify.NIRSmat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{6}.spm.tools.nirs10.model_specify.wls_bglm_specify.NIRSmat(1).sname = 'NIRScoreg: NIRS.mat';
matlabbatch{6}.spm.tools.nirs10.model_specify.wls_bglm_specify.NIRSmat(1).src_exbranch = substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{6}.spm.tools.nirs10.model_specify.wls_bglm_specify.NIRSmat(1).src_output = substruct('.','NIRSmat');
matlabbatch{6}.spm.tools.nirs10.model_specify.wls_bglm_specify.dir1 = 'Stat';
matlabbatch{6}.spm.tools.nirs10.model_specify.wls_bglm_specify.subj.input_onsets = {''};
matlabbatch{6}.spm.tools.nirs10.model_specify.wls_bglm_specify.subj.multi_reg = {''};
matlabbatch{6}.spm.tools.nirs10.model_specify.wls_bglm_specify.units = 1;
matlabbatch{6}.spm.tools.nirs10.model_specify.wls_bglm_specify.time_res = 1;
matlabbatch{6}.spm.tools.nirs10.model_specify.wls_bglm_specify.derivs = [0 0];
matlabbatch{6}.spm.tools.nirs10.model_specify.wls_bglm_specify.bases.hrf.derivs = [0 0];
matlabbatch{6}.spm.tools.nirs10.model_specify.wls_bglm_specify.volt = 1;
matlabbatch{6}.spm.tools.nirs10.model_specify.wls_bglm_specify.GLM_include_cardiac = 0;
matlabbatch{6}.spm.tools.nirs10.model_specify.wls_bglm_specify.GLM_include_Mayer = 0;
matlabbatch{6}.spm.tools.nirs10.model_specify.wls_bglm_specify.NIRSchannelsConfound.NoNIRSconfounds = struct([]);
matlabbatch{6}.spm.tools.nirs10.model_specify.wls_bglm_specify.GenerateHbT = 1;
matlabbatch{6}.spm.tools.nirs10.model_specify.wls_bglm_specify.flag_window = 0;
matlabbatch{6}.spm.tools.nirs10.model_specify.wls_bglm_specify.channel_pca = 0;
matlabbatch{6}.spm.tools.nirs10.model_specify.wls_bglm_specify.hpf_butter.hpf_butter_Off = struct([]);
matlabbatch{6}.spm.tools.nirs10.model_specify.wls_bglm_specify.generate_trRV = 1;
matlabbatch{6}.spm.tools.nirs10.model_specify.wls_bglm_specify.filter_design_matrix = 1;
matlabbatch{6}.spm.tools.nirs10.model_specify.wls_bglm_specify.wls_or_bglm.NIRS_SPM.nirs_noise = 0;
matlabbatch{6}.spm.tools.nirs10.model_specify.wls_bglm_specify.wls_or_bglm.NIRS_SPM.nirs_hpf.hpf_wavelet.hpf_wavelet_iter = 4;
matlabbatch{6}.spm.tools.nirs10.model_specify.wls_bglm_specify.wls_or_bglm.NIRS_SPM.nirs_lpf.lpf_gauss.fwhm1 = 1.5;
matlabbatch{6}.spm.tools.nirs10.model_specify.wls_bglm_specify.LiomDeleteLarge = 0;
matlabbatch{7}.spm.tools.nirs10.model_estimate.wls_bglm_estimate.NIRSmat(1) = cfg_dep;
matlabbatch{7}.spm.tools.nirs10.model_estimate.wls_bglm_estimate.NIRSmat(1).tname = 'NIRS.mat';
matlabbatch{7}.spm.tools.nirs10.model_estimate.wls_bglm_estimate.NIRSmat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{7}.spm.tools.nirs10.model_estimate.wls_bglm_estimate.NIRSmat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{7}.spm.tools.nirs10.model_estimate.wls_bglm_estimate.NIRSmat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{7}.spm.tools.nirs10.model_estimate.wls_bglm_estimate.NIRSmat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{7}.spm.tools.nirs10.model_estimate.wls_bglm_estimate.NIRSmat(1).sname = 'LIOM GLM Specification: NIRS.mat';
matlabbatch{7}.spm.tools.nirs10.model_estimate.wls_bglm_estimate.NIRSmat(1).src_exbranch = substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{7}.spm.tools.nirs10.model_estimate.wls_bglm_estimate.NIRSmat(1).src_output = substruct('.','NIRSmat');
matlabbatch{7}.spm.tools.nirs10.model_estimate.wls_bglm_estimate.NIRS_SPM_which_GLM = 1;
matlabbatch{8}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1) = cfg_dep;
matlabbatch{8}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).tname = 'NIRS.mat';
matlabbatch{8}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{8}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{8}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{8}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{8}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).sname = 'LIOM GLM Estimation: NIRS.mat';
matlabbatch{8}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).src_exbranch = substruct('.','val', '{}',{7}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{8}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).src_output = substruct('.','NIRSmat');
matlabbatch{8}.spm.tools.nirs10.model_estimate.liom_contrast.NewDirCopyNIRS.CreateNIRSCopy_false = struct([]);
matlabbatch{8}.spm.tools.nirs10.model_estimate.liom_contrast.view = [5 3 4];
matlabbatch{8}.spm.tools.nirs10.model_estimate.liom_contrast.consess = {};
matlabbatch{8}.spm.tools.nirs10.model_estimate.liom_contrast.spatial_LPF.spatial_LPF_Off = struct([]);
matlabbatch{8}.spm.tools.nirs10.model_estimate.liom_contrast.ProcessContrastsBySession = 1;
matlabbatch{8}.spm.tools.nirs10.model_estimate.liom_contrast.GenerateInverted = 1;
matlabbatch{8}.spm.tools.nirs10.model_estimate.liom_contrast.GroupColorbars = 0;
matlabbatch{8}.spm.tools.nirs10.model_estimate.liom_contrast.contrast_p_value = 0.05;
matlabbatch{8}.spm.tools.nirs10.model_estimate.liom_contrast.contrast_figures = 3;
matlabbatch{8}.spm.tools.nirs10.model_estimate.liom_contrast.override_colorbar.colorbar_default = struct([]);
matlabbatch{8}.spm.tools.nirs10.model_estimate.liom_contrast.figures_visible = 0;
matlabbatch{8}.spm.tools.nirs10.model_estimate.liom_contrast.GroupFiguresIntoSubplots = 1;
matlabbatch{8}.spm.tools.nirs10.model_estimate.liom_contrast.output_unc = 0;
matlabbatch{8}.spm.tools.nirs10.model_estimate.liom_contrast.SmallFigures = 1;
matlabbatch{8}.spm.tools.nirs10.model_estimate.liom_contrast.write_neg_pos = 0;
matlabbatch{8}.spm.tools.nirs10.model_estimate.liom_contrast.TopoData = {''};
matlabbatch{8}.spm.tools.nirs10.model_estimate.liom_contrast.GroupMultiSession = 0;
matlabbatch{9}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1) = cfg_dep;
matlabbatch{9}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).tname = 'NIRS.mat';
matlabbatch{9}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{9}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{9}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{9}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{9}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).sname = 'Liom Contrast Calculations: NIRS.mat';
matlabbatch{9}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).src_exbranch = substruct('.','val', '{}',{8}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{9}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).src_output = substruct('.','NIRSmat');
matlabbatch{9}.spm.tools.nirs10.model_estimate.liom_group.FFX_or_RFX = 0;
matlabbatch{9}.spm.tools.nirs10.model_estimate.liom_group.contrast_figures = 3;
matlabbatch{9}.spm.tools.nirs10.model_estimate.liom_group.contrast_p_value = 0.05;
matlabbatch{9}.spm.tools.nirs10.model_estimate.liom_group.GenerateInverted = 1;
matlabbatch{9}.spm.tools.nirs10.model_estimate.liom_group.GroupColorbars = 0;
matlabbatch{9}.spm.tools.nirs10.model_estimate.liom_group.override_colorbar.colorbar_default = struct([]);
matlabbatch{9}.spm.tools.nirs10.model_estimate.liom_group.figures_visible = 0;
matlabbatch{9}.spm.tools.nirs10.model_estimate.liom_group.GroupFiguresIntoSubplots = 1;
matlabbatch{9}.spm.tools.nirs10.model_estimate.liom_group.output_unc = 0;
matlabbatch{9}.spm.tools.nirs10.model_estimate.liom_group.SmallFigures = 1;
matlabbatch{9}.spm.tools.nirs10.model_estimate.liom_group.write_neg_pos = 0;
matlabbatch{9}.spm.tools.nirs10.model_estimate.liom_group.group_session_to_average = 1;
matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1) = cfg_dep;
matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).tname = 'NIRS.mat';
matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).sname = 'Liom Group Model Estimation: NIRS.mat';
matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).src_exbranch = substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).src_output = substruct('.','NIRSmat');
matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_group.FFX_or_RFX = 0;
matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_group.contrast_figures = 3;
matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_group.contrast_p_value = 0.01;
matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_group.GenerateInverted = 1;
matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_group.GroupColorbars = 0;
matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_group.override_colorbar.colorbar_default = struct([]);
matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_group.figures_visible = 0;
matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_group.GroupFiguresIntoSubplots = 1;
matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_group.output_unc = 0;
matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_group.SmallFigures = 1;
matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_group.write_neg_pos = 0;
matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_group.group_session_to_average = 1;
matlabbatch{11}.spm.tools.nirs10.model_estimate.extract_map_data.NIRSmat(1) = cfg_dep;
matlabbatch{11}.spm.tools.nirs10.model_estimate.extract_map_data.NIRSmat(1).tname = 'NIRS.mat';
matlabbatch{11}.spm.tools.nirs10.model_estimate.extract_map_data.NIRSmat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{11}.spm.tools.nirs10.model_estimate.extract_map_data.NIRSmat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{11}.spm.tools.nirs10.model_estimate.extract_map_data.NIRSmat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{11}.spm.tools.nirs10.model_estimate.extract_map_data.NIRSmat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{11}.spm.tools.nirs10.model_estimate.extract_map_data.NIRSmat(1).sname = 'Liom Group Model Estimation: NIRS.mat';
matlabbatch{11}.spm.tools.nirs10.model_estimate.extract_map_data.NIRSmat(1).src_exbranch = substruct('.','val', '{}',{10}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{11}.spm.tools.nirs10.model_estimate.extract_map_data.NIRSmat(1).src_output = substruct('.','NIRSmat');
matlabbatch{11}.spm.tools.nirs10.model_estimate.extract_map_data.view = 5;
matlabbatch{11}.spm.tools.nirs10.model_estimate.extract_map_data.extract_contrast = 1;
matlabbatch{11}.spm.tools.nirs10.model_estimate.extract_map_data.extract_select_mode.extract_auto_mode.extract_select_auto_mode.extract_max_HbR = struct([]);
matlabbatch{11}.spm.tools.nirs10.model_estimate.extract_map_data.extract_average_mode.extract_threshold.extract_threshold_val = 3.5;
matlabbatch{11}.spm.tools.nirs10.model_estimate.extract_map_data.extract_average_mode.extract_threshold.extract_radius_val = 3;
matlabbatch{11}.spm.tools.nirs10.model_estimate.extract_map_data.extract_struct_name = 'ED';