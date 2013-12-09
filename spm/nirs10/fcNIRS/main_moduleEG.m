function main_moduleEG(OP)
subject = OP.subjects;
mux = OP.mux;
subj = OP.subj;
path0 = OP.path0;
path1 = OP.path1;
pathbatch = OP.pathbatch;
for s1=1:length(subject)
    try
    n0 = subject(s1);
    cs = subj{n0};
    dp = fullfile(path0,cs);%data path
    cp = fullfile(path1,cs);%current path
    clear dir_T1 dir_EEG dir_MTG dir_DATA
    %0) find the subject directory
    [filesRec,DirsRec] = spm_select('FPListRec',dp,'.*');
    [filesRec_cp,DirsRec_cp] = spm_select('FPListRec',cp,'.*');
    %0.1) directory match
    for i0=1:length(DirsRec)
        %T1
        if strfind(DirsRec{i0},'T1')
            dir_T1 = DirsRec{i0};
        elseif strfind(DirsRec{i0},'eeg')
            dir_EEG = DirsRec{i0};
        elseif strfind(DirsRec{i0},'mtg')
            dir_MTG = DirsRec{i0};
        elseif strfind(DirsRec{i0},'dataBOXY')
            dir_DATA = DirsRec{i0};
        else
            continue;
        end
    end
    
    for i0 = 1 : length(DirsRec_cp)
        if strfind(DirsRec_cp{i0},'T1')
            dir_T1 = DirsRec_cp{i0};
        end
    end
    
    %1) find BOXYfiles
    [filesRec_dataBOXY,DirsRec_dataBOXY] = spm_select('FPListRec',dir_DATA,'.*');
    %2) find EEG files
%     [filesRec_eeg,DirsRec_eeg] = spm_select('FPListRec',dir_EEG,'.*');
    %3) find mtg files
    [filesRec_mtg,DirsRec_mtg] = spm_select('FPListRec',dir_MTG,'.prj');
    %4) find T1 file(s)
%     [filesRec_T1,DirsRec_T1] = spm_select('FPList',dir_T1,'epi*');
    
    if OP.first_session
        filesRec_dataBOXY = {filesRec_dataBOXY{1}};
%         filesRec_eeg = {filesRec_eeg{1}};
    end
    
    %To search for the correct T1 image
    found = [];
    
%     for f0=1:size(filesRec_T1,1)
%         file_T1 = filesRec_T1(f0,:);
%         [dir1 fil1 ext1] = fileparts(file_T1);
%         if strcmp(fil1(1:3),'epi')
%             %should be a good c1 nii file
%             found = file_T1;
%             break;
%         end
%     end
    found_T1 = [];
    for f0=1:size(found,1)
        file_T1 = found(f0,:);
        [dir1 fil1 ext1] = fileparts(file_T1);
        if strcmp(deblank(ext1),'.nii')
            %should be a good c1 nii file
            found_T1 = deblank(file_T1);
            break;
        end
    end
    
    if ~isempty(found_T1)
        filesRec_T1 = found_T1;
    else
        %No render image found
%         disp(['ERROR! No T1 image found under the directory' dir_T1]);
    end
    
    if ~isempty(found_T1)
        [dir_T1_found file_T1_found ext_T1_found] = fileparts(found);
        filesRec_wT1 = fullfile(dir_T1_found,['w' file_T1_found ext_T1_found]);
        filesRec_c1 = fullfile(dir_T1_found,['c1' file_T1_found ext_T1_found]);
        filesRec_wc1 = fullfile(dir_T1_found,['wc1' file_T1_found ext_T1_found]);
    end
    cp_full = fullfile(cp);
    clear matlabbatch
    matlabbatch{1}.spm.tools.nirs10.readNIRS.boxy_manual1.subj2.Apath = {cp_full};
    matlabbatch{1}.spm.tools.nirs10.readNIRS.boxy_manual1.subj2.fnames = filesRec_dataBOXY;
    matlabbatch{1}.spm.tools.nirs10.readNIRS.boxy_manual1.subj2.prjFile = filesRec_mtg;
    matlabbatch{1}.spm.tools.nirs10.readNIRS.boxy_manual1.subj2.age1 = 25;
    matlabbatch{1}.spm.tools.nirs10.readNIRS.boxy_manual1.subj2.raw_onset_files2 = filesRec_eeg;
    matlabbatch{1}.spm.tools.nirs10.readNIRS.boxy_manual1.subj2.anatT1 =  {[filesRec_T1 ',1']};
    matlabbatch{1}.spm.tools.nirs10.readNIRS.boxy_manual1.force_redo = 0;
    matlabbatch{1}.spm.tools.nirs10.readNIRS.boxy_manual1.config_path2.T1_path = 'T1';
    matlabbatch{1}.spm.tools.nirs10.readNIRS.boxy_manual1.config_path2.output_path = OP.main_dir;
    matlabbatch{1}.spm.tools.nirs10.readNIRS.boxy_manual1.cf1.Lambda = [830 690];
    matlabbatch{1}.spm.tools.nirs10.readNIRS.boxy_manual1.cf1.freq = 19.5312;
    matlabbatch{1}.spm.tools.nirs10.readNIRS.boxy_manual1.cf1.distmin = 1;
    matlabbatch{1}.spm.tools.nirs10.readNIRS.boxy_manual1.cf1.distmax = 6;
    matlabbatch{1}.spm.tools.nirs10.readNIRS.boxy_manual1.cf1.save_bin1 = true;
    matlabbatch{1}.spm.tools.nirs10.readNIRS.boxy_manual1.cf1.sizebloc = 1024;
    matlabbatch{1}.spm.tools.nirs10.readNIRS.boxy_manual1.cf1.nb_Mux = mux{n0};
    matlabbatch{1}.spm.tools.nirs10.readNIRS.boxy_manual1.cf1.MaxSources = 64;
    matlabbatch{1}.spm.tools.nirs10.readNIRS.boxy_manual1.cf1.nb_Det = 16;
    matlabbatch{1}.spm.tools.nirs10.readNIRS.boxy_manual1.cf1.MaxElectrodes = 19;
    matlabbatch{1}.spm.tools.nirs10.readNIRS.boxy_manual1.cf1.use10_10system = true;
    matlabbatch{1}.spm.tools.nirs10.readNIRS.boxy_manual1.cf1.resample = 1;
    matlabbatch{1}.spm.tools.nirs10.readNIRS.boxy_manual1.cf1.sample_LPF.sample_LPF_off = struct([]);
    matlabbatch{2}.spm.tools.nirs10.preprocANAT.MCsegment1.NIRSmat(1) = cfg_dep;
    matlabbatch{2}.spm.tools.nirs10.preprocANAT.MCsegment1.NIRSmat(1).tname = 'NIRS.mat (optional)';
    matlabbatch{2}.spm.tools.nirs10.preprocANAT.MCsegment1.NIRSmat(1).tgt_spec{1}(1).name = 'filter';
    matlabbatch{2}.spm.tools.nirs10.preprocANAT.MCsegment1.NIRSmat(1).tgt_spec{1}(1).value = 'mat';
    matlabbatch{2}.spm.tools.nirs10.preprocANAT.MCsegment1.NIRSmat(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{2}.spm.tools.nirs10.preprocANAT.MCsegment1.NIRSmat(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{2}.spm.tools.nirs10.preprocANAT.MCsegment1.NIRSmat(1).sname = 'ReadBoxyManual: NIRS.mat';
    matlabbatch{2}.spm.tools.nirs10.preprocANAT.MCsegment1.NIRSmat(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{2}.spm.tools.nirs10.preprocANAT.MCsegment1.NIRSmat(1).src_output = substruct('.','NIRSmat');
    matlabbatch{2}.spm.tools.nirs10.preprocANAT.MCsegment1.force_redo = 0;
    matlabbatch{2}.spm.tools.nirs10.preprocANAT.MCsegment1.NIRSmatCopyChoice.NIRSmatOverwrite = struct([]);
    matlabbatch{2}.spm.tools.nirs10.preprocANAT.MCsegment1.force_reprocess = 0;
    matlabbatch{2}.spm.tools.nirs10.preprocANAT.MCsegment1.image_in = {''};
    matlabbatch{2}.spm.tools.nirs10.preprocANAT.MCsegment1.output_autonaming = 0;
    matlabbatch{2}.spm.tools.nirs10.preprocANAT.MCsegment1.output_prefix = 'Only edit if you chose ''No'' to ''Automatic output naming''';
    matlabbatch{2}.spm.tools.nirs10.preprocANAT.MCsegment1.skn.sorting_method = 1;
    matlabbatch{2}.spm.tools.nirs10.preprocANAT.MCsegment1.skl.sorting_method = 2;
    matlabbatch{2}.spm.tools.nirs10.preprocANAT.MCsegment1.csf.sorting_method = 0;
    matlabbatch{2}.spm.tools.nirs10.preprocANAT.MCsegment1.grm.sorting_method = 0;
    matlabbatch{2}.spm.tools.nirs10.preprocANAT.MCsegment1.wtm.sorting_method = 0;
    matlabbatch{2}.spm.tools.nirs10.preprocANAT.MCsegment1.vbm_seg = 0;
    matlabbatch{2}.spm.tools.nirs10.preprocANAT.MCsegment1.thresh_as = 0.6;
    matlabbatch{2}.spm.tools.nirs10.preprocANAT.MCsegment1.head_shadow.thresh_hs = 0.6;
    matlabbatch{2}.spm.tools.nirs10.preprocANAT.MCsegment1.head_shadow.se_size_hs = 2;
    matlabbatch{2}.spm.tools.nirs10.preprocANAT.MCsegment1.rebel_surrounding = 3;
    matlabbatch{2}.spm.tools.nirs10.preprocANAT.MCsegment1.rebel_thresh_hs = 0.3;
    matlabbatch{2}.spm.tools.nirs10.preprocANAT.MCsegment1.process_image.se_size_pi = 2;
    matlabbatch{2}.spm.tools.nirs10.preprocANAT.MCsegment1.process_image.gaussfilt_size = 7;
    matlabbatch{2}.spm.tools.nirs10.preprocANAT.MCsegment1.process_image.gaussfilt_sdev = 4;
    matlabbatch{3}.spm.tools.nirs10.readOnsets.AnalyzerOnsets.NIRSmat(1) = cfg_dep;
    matlabbatch{3}.spm.tools.nirs10.readOnsets.AnalyzerOnsets.NIRSmat(1).tname = 'NIRS.mat (optional)';
    matlabbatch{3}.spm.tools.nirs10.readOnsets.AnalyzerOnsets.NIRSmat(1).tgt_spec{1}(1).name = 'filter';
    matlabbatch{3}.spm.tools.nirs10.readOnsets.AnalyzerOnsets.NIRSmat(1).tgt_spec{1}(1).value = 'mat';
    matlabbatch{3}.spm.tools.nirs10.readOnsets.AnalyzerOnsets.NIRSmat(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{3}.spm.tools.nirs10.readOnsets.AnalyzerOnsets.NIRSmat(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{3}.spm.tools.nirs10.readOnsets.AnalyzerOnsets.NIRSmat(1).sname = 'MC Segmentation: NIRS.mat';
    matlabbatch{3}.spm.tools.nirs10.readOnsets.AnalyzerOnsets.NIRSmat(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{3}.spm.tools.nirs10.readOnsets.AnalyzerOnsets.NIRSmat(1).src_output = substruct('.','NIRSmat');
    matlabbatch{3}.spm.tools.nirs10.readOnsets.AnalyzerOnsets.force_redo = 0;
    matlabbatch{3}.spm.tools.nirs10.readOnsets.AnalyzerOnsets.NIRSmatCopyChoice.NIRSmatOverwrite = struct([]);
    matlabbatch{3}.spm.tools.nirs10.readOnsets.AnalyzerOnsets.raw_onset_files = {''};
    matlabbatch{3}.spm.tools.nirs10.readOnsets.AnalyzerOnsets.onset_to_keep = '';
    matlabbatch{3}.spm.tools.nirs10.readOnsets.AnalyzerOnsets.freq_NIRS1 = [];
    matlabbatch{3}.spm.tools.nirs10.readOnsets.AnalyzerOnsets.dp_NIRS1 = [];
    matlabbatch{3}.spm.tools.nirs10.readOnsets.AnalyzerOnsets.cardiac_repair.cardiac_repair_on.avg_number = 5;
    matlabbatch{3}.spm.tools.nirs10.readOnsets.AnalyzerOnsets.cardiac_repair.cardiac_repair_on.gap_def = 1.8;
    matlabbatch{4}.spm.tools.nirs10.preprocessNIRS.remove_chn_stdev.NIRSmat(1) = cfg_dep;
    matlabbatch{4}.spm.tools.nirs10.preprocessNIRS.remove_chn_stdev.NIRSmat(1).tname = 'NIRS.mat';
    matlabbatch{4}.spm.tools.nirs10.preprocessNIRS.remove_chn_stdev.NIRSmat(1).tgt_spec{1}(1).name = 'filter';
    matlabbatch{4}.spm.tools.nirs10.preprocessNIRS.remove_chn_stdev.NIRSmat(1).tgt_spec{1}(1).value = 'mat';
    matlabbatch{4}.spm.tools.nirs10.preprocessNIRS.remove_chn_stdev.NIRSmat(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{4}.spm.tools.nirs10.preprocessNIRS.remove_chn_stdev.NIRSmat(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{4}.spm.tools.nirs10.preprocessNIRS.remove_chn_stdev.NIRSmat(1).sname = 'Read NIRS onsets: NIRS.mat';
    matlabbatch{4}.spm.tools.nirs10.preprocessNIRS.remove_chn_stdev.NIRSmat(1).src_exbranch = substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{4}.spm.tools.nirs10.preprocessNIRS.remove_chn_stdev.NIRSmat(1).src_output = substruct('.','NIRSmat');
    matlabbatch{4}.spm.tools.nirs10.preprocessNIRS.remove_chn_stdev.force_redo = 0;
    matlabbatch{4}.spm.tools.nirs10.preprocessNIRS.remove_chn_stdev.NIRSmatCopyChoice.NIRSmatOverwrite = struct([]);
    matlabbatch{4}.spm.tools.nirs10.preprocessNIRS.remove_chn_stdev.threshold_stdev = 0.2;
    matlabbatch{4}.spm.tools.nirs10.preprocessNIRS.remove_chn_stdev.window_stdev = 5;
    matlabbatch{5}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.NIRSmat(1) = cfg_dep;
    matlabbatch{5}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.NIRSmat(1).tname = 'NIRS.mat';
    matlabbatch{5}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.NIRSmat(1).tgt_spec{1}(1).name = 'filter';
    matlabbatch{5}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.NIRSmat(1).tgt_spec{1}(1).value = 'mat';
    matlabbatch{5}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.NIRSmat(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{5}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.NIRSmat(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{5}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.NIRSmat(1).sname = 'Remove noisy channels (stdev): NIRS.mat';
    matlabbatch{5}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.NIRSmat(1).src_exbranch = substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{5}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.NIRSmat(1).src_output = substruct('.','NIRSmat');
    matlabbatch{5}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.force_redo = 0;
    matlabbatch{5}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.NIRSmatCopyChoice.NIRSmatOverwrite = struct([]);
    matlabbatch{5}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.Normalize_OD = 0;
    matlabbatch{5}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.add_or_mult = 0;
    matlabbatch{5}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.baseline_duration = 2;
    matlabbatch{5}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.normalization_type = 2;
    matlabbatch{5}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.Analyzer_sf = 1;
    matlabbatch{5}.spm.tools.nirs10.preprocessNIRS.normalize_baseline.nirs_filling_jumps.nirs_fill_jumps_off = struct([]);
    matlabbatch{6}.spm.tools.nirs10.preprocessNIRS.ODtoHbOHbR.NIRSmat(1) = cfg_dep;
    matlabbatch{6}.spm.tools.nirs10.preprocessNIRS.ODtoHbOHbR.NIRSmat(1).tname = 'NIRS.mat';
    matlabbatch{6}.spm.tools.nirs10.preprocessNIRS.ODtoHbOHbR.NIRSmat(1).tgt_spec{1}(1).name = 'filter';
    matlabbatch{6}.spm.tools.nirs10.preprocessNIRS.ODtoHbOHbR.NIRSmat(1).tgt_spec{1}(1).value = 'mat';
    matlabbatch{6}.spm.tools.nirs10.preprocessNIRS.ODtoHbOHbR.NIRSmat(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{6}.spm.tools.nirs10.preprocessNIRS.ODtoHbOHbR.NIRSmat(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{6}.spm.tools.nirs10.preprocessNIRS.ODtoHbOHbR.NIRSmat(1).sname = 'Normalize Baseline: NIRS.mat';
    matlabbatch{6}.spm.tools.nirs10.preprocessNIRS.ODtoHbOHbR.NIRSmat(1).src_exbranch = substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{6}.spm.tools.nirs10.preprocessNIRS.ODtoHbOHbR.NIRSmat(1).src_output = substruct('.','NIRSmat');
    matlabbatch{6}.spm.tools.nirs10.preprocessNIRS.ODtoHbOHbR.force_redo = 0;
    matlabbatch{6}.spm.tools.nirs10.preprocessNIRS.ODtoHbOHbR.NIRSmatCopyChoice.NIRSmatOverwrite = struct([]);
    matlabbatch{6}.spm.tools.nirs10.preprocessNIRS.ODtoHbOHbR.DPF.DPFlit = struct([]);
    matlabbatch{6}.spm.tools.nirs10.preprocessNIRS.ODtoHbOHbR.PVF.PVFval = [50 50];
    matlabbatch{6}.spm.tools.nirs10.preprocessNIRS.ODtoHbOHbR.nirs_filling_jumps.nirs_fill_jumps_off = struct([]);
    if OP.do_coreg
        matlabbatch{7}.spm.tools.nirs10.coregNIRS.coregnew1.NIRSmat(1) = cfg_dep;
        matlabbatch{7}.spm.tools.nirs10.coregNIRS.coregnew1.NIRSmat(1).tname = 'NIRS.mat';
        matlabbatch{7}.spm.tools.nirs10.coregNIRS.coregnew1.NIRSmat(1).tgt_spec{1}(1).name = 'filter';
        matlabbatch{7}.spm.tools.nirs10.coregNIRS.coregnew1.NIRSmat(1).tgt_spec{1}(1).value = 'mat';
        matlabbatch{7}.spm.tools.nirs10.coregNIRS.coregnew1.NIRSmat(1).tgt_spec{1}(2).name = 'strtype';
        matlabbatch{7}.spm.tools.nirs10.coregNIRS.coregnew1.NIRSmat(1).tgt_spec{1}(2).value = 'e';
        matlabbatch{7}.spm.tools.nirs10.coregNIRS.coregnew1.NIRSmat(1).sname = 'Convert OD to HbO/HbR : NIRS.mat';
        matlabbatch{7}.spm.tools.nirs10.coregNIRS.coregnew1.NIRSmat(1).src_exbranch = substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
        matlabbatch{7}.spm.tools.nirs10.coregNIRS.coregnew1.NIRSmat(1).src_output = substruct('.','NIRSmat');
        matlabbatch{7}.spm.tools.nirs10.coregNIRS.coregnew1.force_redo = OP.force_coreg;
        matlabbatch{7}.spm.tools.nirs10.coregNIRS.coregnew1.ForceReprocess = 0;
        matlabbatch{7}.spm.tools.nirs10.coregNIRS.coregnew1.NIRSmatCopyChoice.NIRSmatCopy.NewNIRSdir = OP.coreg_dir;
        matlabbatch{7}.spm.tools.nirs10.coregNIRS.coregnew1.fiducial_MNI_choice = OP.subj_fiducials{n0};
        matlabbatch{7}.spm.tools.nirs10.coregNIRS.coregnew1.radius_channel = OP.coreg_size;
        if OP.subj_fiducials{n0}
            Fcoord = OP.Fcoord{n0};
            matlabbatch{7}.spm.tools.nirs10.coregNIRS.coregnew1.nasion_wMNI = Fcoord(1,:); %nasion_MNI{n0};
            matlabbatch{7}.spm.tools.nirs10.coregNIRS.coregnew1.AL_wMNI = Fcoord(2,:); %AL_MNI{n0};
            matlabbatch{7}.spm.tools.nirs10.coregNIRS.coregnew1.AR_wMNI = Fcoord(3,:); %AR_MNI{n0};
        else
            matlabbatch{7}.spm.tools.nirs10.coregNIRS.coregnew1.nasion_wMNI = [0 84 -48];
            matlabbatch{7}.spm.tools.nirs10.coregNIRS.coregnew1.AL_wMNI = [-83 -19 -38];
            matlabbatch{7}.spm.tools.nirs10.coregNIRS.coregnew1.AR_wMNI = [83 -19 -38];
        end
        matlabbatch{7}.spm.tools.nirs10.coregNIRS.coregnew1.cortex_projection_method.project_Korean = struct([]);
        if OP.project_subject{n0}
            matlabbatch{7}.spm.tools.nirs10.coregNIRS.coregnew1.render_choice.render_subject = struct([]);
        else
            matlabbatch{7}.spm.tools.nirs10.coregNIRS.coregnew1.render_choice.render_template = struct([]);
        end
        matlabbatch{7}.spm.tools.nirs10.coregNIRS.coregnew1.OutputSkinFigs = 1;
    end
    if ~OP.onlyCoregandPreprocess
        matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.NIRSmat(1) = cfg_dep;
        matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.NIRSmat(1).tname = 'NIRS.mat';
        matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.NIRSmat(1).tgt_spec{1}(1).name = 'filter';
        matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.NIRSmat(1).tgt_spec{1}(1).value = 'mat';
        matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.NIRSmat(1).tgt_spec{1}(2).name = 'strtype';
        matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.NIRSmat(1).tgt_spec{1}(2).value = 'e';
        matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.NIRSmat(1).sname = 'NIRScoreg (new): NIRS.mat';
        matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.NIRSmat(1).src_exbranch = substruct('.','val', '{}',{7}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
        matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.NIRSmat(1).src_output = substruct('.','NIRSmat');
        matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.force_redo = 0;
        if strcmp(OP.stat_dir,'Stat')
            matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.NIRSmatCopyChoice.NIRSmatOverwrite = struct([]);
        else
            matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.NIRSmatCopyChoice.NIRSmatCopy.NewNIRSdir = OP.stat_dir;
        end
        matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.sessions = '';
        matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.subj.input_onsets = {''};
        matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.subj.multi_reg = {''};
        matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.units = 1;
        matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.time_res = 1;
        if OP.test_hrf
            matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.derivs = [1 1];
        else
            matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.derivs = [0 0];
        end
        if OP.test_hrf
            matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.bases.hrf.derivs = [1 1];
        else
            matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.bases.hrf.derivs = [0 0];
        end
        matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.volt = OP.volt{n0};
        matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.GLM_include_cardiac = OP.includeHR;
        matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.GLM_include_Mayer = 0;
        matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.vasomotion_choice.no_vasomotion = struct([]);
        matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.NIRSchannelsConfound.NoNIRSconfounds = struct([]);
        matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.GenerateHbT = 1;
        matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.flag_window = 1;
        matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.channel_pca = OP.channel_pca{n0};
        matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.NumPCAComponents = OP.nComponents{n0};
        matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.hpf_butter.hpf_butter_On.hpf_butter_freq = OP.HPFfreq;
        matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.hpf_butter.hpf_butter_On.hpf_butter_order = OP.HPForder;
        matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.generate_trRV = 1;
        matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.TrRVRVexact = 1;
        matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.wls_or_bglm.NIRS_SPM.nirs_noise = 0;
        matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.wls_or_bglm.NIRS_SPM.nirs_hpf.hpf_none = struct([]);
        if OP.LPFhrf
            matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.wls_or_bglm.NIRS_SPM.nirs_lpf.lpf_hrf = struct([]);
        else
            matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.wls_or_bglm.NIRS_SPM.nirs_lpf.lpf_gauss.fwhm1 = 1.5;
        end
        if OP.downsize
            matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.target_sampling_rate.specified_sampling_rate.sampling_rate = 1;
        else
            matlabbatch{8}.spm.tools.nirs10.model_specify.wls_bglm_specify.target_sampling_rate.raw_data_sampling_rate = struct([]);
        end
        matlabbatch{9}.spm.tools.nirs10.model_estimate.wls_bglm_estimate.NIRSmat(1) = cfg_dep;
        matlabbatch{9}.spm.tools.nirs10.model_estimate.wls_bglm_estimate.NIRSmat(1).tname = 'NIRS.mat';
        matlabbatch{9}.spm.tools.nirs10.model_estimate.wls_bglm_estimate.NIRSmat(1).tgt_spec{1}(1).name = 'filter';
        matlabbatch{9}.spm.tools.nirs10.model_estimate.wls_bglm_estimate.NIRSmat(1).tgt_spec{1}(1).value = 'mat';
        matlabbatch{9}.spm.tools.nirs10.model_estimate.wls_bglm_estimate.NIRSmat(1).tgt_spec{1}(2).name = 'strtype';
        matlabbatch{9}.spm.tools.nirs10.model_estimate.wls_bglm_estimate.NIRSmat(1).tgt_spec{1}(2).value = 'e';
        matlabbatch{9}.spm.tools.nirs10.model_estimate.wls_bglm_estimate.NIRSmat(1).sname = 'LIOM GLM Specification: NIRS.mat';
        matlabbatch{9}.spm.tools.nirs10.model_estimate.wls_bglm_estimate.NIRSmat(1).src_exbranch = substruct('.','val', '{}',{8}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
        matlabbatch{9}.spm.tools.nirs10.model_estimate.wls_bglm_estimate.NIRSmat(1).src_output = substruct('.','NIRSmat');
        matlabbatch{9}.spm.tools.nirs10.model_estimate.wls_bglm_estimate.force_redo = 0;
        matlabbatch{9}.spm.tools.nirs10.model_estimate.wls_bglm_estimate.NIRSmatCopyChoice.NIRSmatOverwrite = struct([]);
        matlabbatch{9}.spm.tools.nirs10.model_estimate.wls_bglm_estimate.NIRS_SPM_which_GLM = 1;
        matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1) = cfg_dep;
        matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).tname = 'NIRS.mat';
        matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).tgt_spec{1}(1).name = 'filter';
        matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).tgt_spec{1}(1).value = 'mat';
        matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).tgt_spec{1}(2).name = 'strtype';
        matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).tgt_spec{1}(2).value = 'e';
        matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).sname = 'LIOM GLM Estimation: NIRS.mat';
        matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).src_exbranch = substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
        matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).src_output = substruct('.','NIRSmat');
        matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_contrast.force_redo = 0;
        matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmatCopyChoice.NIRSmatOverwrite = struct([]);
        matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_contrast.ContrastChoice.automated_contrasts.NonlinearEpilepsyOn = 0;
        matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_contrast.Sessions = '';
        matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_contrast.view = OP.views{n0};
        matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_contrast.TopoData = {''};
        matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_contrast.StatMethod = 1;
        matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_contrast.UseCorrelRes = 1;
        matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_contrast.spatial_LPF.spatial_LPF_Off = struct([]);
        matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_contrast.contrast_p_value = 0.05;
        matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_contrast.display_options.GenerateInverted = 1;
        matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_contrast.display_options.GroupColorbars = 0;
        matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_contrast.display_options.contrast_figures = 3;
        matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_contrast.display_options.override_colorbar.colorbar_default = struct([]);
        matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_contrast.display_options.figures_visible = 0;
        matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_contrast.display_options.GroupFiguresIntoSubplots = 1;
        matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_contrast.display_options.output_unc = 1;
        matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_contrast.display_options.SmallFigures = 1;
        matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_contrast.display_options.write_neg_pos = 0;
        matlabbatch{10}.spm.tools.nirs10.model_estimate.liom_contrast.display_options.save_nifti_contrasts = 0;
        matlabbatch{11}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1) = cfg_dep;
        matlabbatch{11}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).tname = 'NIRS.mat';
        matlabbatch{11}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).tgt_spec{1}(1).name = 'filter';
        matlabbatch{11}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).tgt_spec{1}(1).value = 'mat';
        matlabbatch{11}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).tgt_spec{1}(2).name = 'strtype';
        matlabbatch{11}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).tgt_spec{1}(2).value = 'e';
        matlabbatch{11}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).sname = 'Liom Contrast Calculations: NIRS.mat';
        matlabbatch{11}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).src_exbranch = substruct('.','val', '{}',{10}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
        matlabbatch{11}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).src_output = substruct('.','NIRSmat');
        matlabbatch{11}.spm.tools.nirs10.model_estimate.liom_group.force_redo = 1;
        matlabbatch{11}.spm.tools.nirs10.model_estimate.liom_group.NIRSmatCopyChoice.NIRSmatOverwrite = struct([]);
        matlabbatch{11}.spm.tools.nirs10.model_estimate.liom_group.group_dir_name = 'Group';
        matlabbatch{11}.spm.tools.nirs10.model_estimate.liom_group.number_dir_to_remove = 3;
        matlabbatch{11}.spm.tools.nirs10.model_estimate.liom_group.FFX_or_RFX = 1;
        matlabbatch{11}.spm.tools.nirs10.model_estimate.liom_group.ContrastChoice.automated_contrasts = struct([]);
        matlabbatch{11}.spm.tools.nirs10.model_estimate.liom_group.StatMethod = 1;
        matlabbatch{11}.spm.tools.nirs10.model_estimate.liom_group.contrast_p_value = 0.05;
        matlabbatch{11}.spm.tools.nirs10.model_estimate.liom_group.group_session_to_average = 1;
        matlabbatch{11}.spm.tools.nirs10.model_estimate.liom_group.simple_sum = 0;
        matlabbatch{11}.spm.tools.nirs10.model_estimate.liom_group.display_options.GenerateInverted = 1;
        matlabbatch{11}.spm.tools.nirs10.model_estimate.liom_group.display_options.GroupColorbars = 0;
        matlabbatch{11}.spm.tools.nirs10.model_estimate.liom_group.display_options.contrast_figures = 3;
        matlabbatch{11}.spm.tools.nirs10.model_estimate.liom_group.display_options.override_colorbar.colorbar_default = struct([]);
        matlabbatch{11}.spm.tools.nirs10.model_estimate.liom_group.display_options.figures_visible = 0;
        matlabbatch{11}.spm.tools.nirs10.model_estimate.liom_group.display_options.GroupFiguresIntoSubplots = 1;
        matlabbatch{11}.spm.tools.nirs10.model_estimate.liom_group.display_options.output_unc = 1;
        matlabbatch{11}.spm.tools.nirs10.model_estimate.liom_group.display_options.SmallFigures = 1;
        matlabbatch{11}.spm.tools.nirs10.model_estimate.liom_group.display_options.write_neg_pos = 0;
        matlabbatch{11}.spm.tools.nirs10.model_estimate.liom_group.display_options.save_nifti_contrasts = 0;
        
        if OP.intrasubject_average
            matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.NIRSmat(1) = cfg_dep;
            matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.NIRSmat(1) = cfg_dep;
            matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.NIRSmat(1).tname = 'NIRS.mat';
            matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.NIRSmat(1).tgt_spec{1}(1).name = 'filter';
            matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.NIRSmat(1).tgt_spec{1}(1).value = 'mat';
            matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.NIRSmat(1).tgt_spec{1}(2).name = 'strtype';
            matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.NIRSmat(1).tgt_spec{1}(2).value = 'e';
            matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.NIRSmat(1).sname = 'NIRScoreg: NIRS.mat';
            matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.NIRSmat(1).src_exbranch = substruct('.','val', '{}',{7}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
            matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.NIRSmat(1).src_output = substruct('.','NIRSmat');
            matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.force_redo = 0;
            matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.NIRSmatCopyChoice.NIRSmatCopy.NewNIRSdir = 'Avg';
            matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.sessions = '';
            matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.subj.input_onsets = {''};
            matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.subj.multi_reg = {''};
            matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.units = 1;
            matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.time_res = 1;
            matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.GLM_include_cardiac = 0;
            matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.GLM_include_Mayer = 0;
            matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.vasomotion_choice.no_vasomotion = struct([]);
            matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.NIRSchannelsConfound.NoNIRSconfounds = struct([]);
            matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.GenerateHbT = 1;
            matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.channel_pca = 0;
            matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.hpf_butter.hpf_butter_On.hpf_butter_freq = OP.HPFfreq;
            matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.hpf_butter.hpf_butter_On.hpf_butter_order = OP.HPForder;
            %matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.hpf_butter.hpf_butter_Off = struct([]);
            matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.AvgFilters.nirs_hpf.hpf_none = struct([]);
            matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.AvgFilters.nirs_lpf.lpf_gauss.fwhm1 = 1.5;
            matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.averaging_choice.block_averaging.onset_delay = 2;
            matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.averaging_choice.block_averaging.onset_duration = 6;
            matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.baseline_choice.baseline_block_averaging.baseline_offset = 0;
            matlabbatch{12}.spm.tools.nirs10.model_specify.liom_intrasubject_average.baseline_choice.baseline_block_averaging.baseline_duration = 2;
            matlabbatch{13}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1) = cfg_dep;
            matlabbatch{13}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).tname = 'NIRS.mat';
            matlabbatch{13}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).tgt_spec{1}(1).name = 'filter';
            matlabbatch{13}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).tgt_spec{1}(1).value = 'mat';
            matlabbatch{13}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).tgt_spec{1}(2).name = 'strtype';
            matlabbatch{13}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).tgt_spec{1}(2).value = 'e';
            matlabbatch{13}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).sname = 'LIOM Intrasubject Average: NIRS.mat';
            matlabbatch{13}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).src_exbranch = substruct('.','val', '{}',{12}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
            matlabbatch{13}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmat(1).src_output = substruct('.','NIRSmat');
            matlabbatch{13}.spm.tools.nirs10.model_estimate.liom_contrast.force_redo = 0;
            matlabbatch{13}.spm.tools.nirs10.model_estimate.liom_contrast.NIRSmatCopyChoice.NIRSmatOverwrite = struct([]);
            matlabbatch{13}.spm.tools.nirs10.model_estimate.liom_contrast.ContrastChoice.automated_contrasts.NonlinearEpilepsyOn = 0;
            matlabbatch{13}.spm.tools.nirs10.model_estimate.liom_contrast.Sessions = '';
            matlabbatch{13}.spm.tools.nirs10.model_estimate.liom_contrast.view = OP.views{n0};
            matlabbatch{13}.spm.tools.nirs10.model_estimate.liom_contrast.TopoData = {''};
            matlabbatch{13}.spm.tools.nirs10.model_estimate.liom_contrast.StatMethod = 1;
            matlabbatch{13}.spm.tools.nirs10.model_estimate.liom_contrast.UseCorrelRes = 1;
            matlabbatch{13}.spm.tools.nirs10.model_estimate.liom_contrast.spatial_LPF.spatial_LPF_Off = struct([]);
            matlabbatch{13}.spm.tools.nirs10.model_estimate.liom_contrast.contrast_p_value = 0.05;
            matlabbatch{13}.spm.tools.nirs10.model_estimate.liom_contrast.display_options.GenerateInverted = 1;
            matlabbatch{13}.spm.tools.nirs10.model_estimate.liom_contrast.display_options.GroupColorbars = 0;
            matlabbatch{13}.spm.tools.nirs10.model_estimate.liom_contrast.display_options.contrast_figures = 3;
            matlabbatch{13}.spm.tools.nirs10.model_estimate.liom_contrast.display_options.override_colorbar.colorbar_default = struct([]);
            matlabbatch{13}.spm.tools.nirs10.model_estimate.liom_contrast.display_options.figures_visible = 0;
            matlabbatch{13}.spm.tools.nirs10.model_estimate.liom_contrast.display_options.GroupFiguresIntoSubplots = 1;
            matlabbatch{13}.spm.tools.nirs10.model_estimate.liom_contrast.display_options.output_unc = 1;
            matlabbatch{13}.spm.tools.nirs10.model_estimate.liom_contrast.display_options.SmallFigures = 1;
            matlabbatch{13}.spm.tools.nirs10.model_estimate.liom_contrast.display_options.write_neg_pos = 0;
            matlabbatch{13}.spm.tools.nirs10.model_estimate.liom_contrast.display_options.save_nifti_contrasts = 0;
            matlabbatch{14}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1) = cfg_dep;
            matlabbatch{14}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).tname = 'NIRS.mat';
            matlabbatch{14}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).tgt_spec{1}(1).name = 'filter';
            matlabbatch{14}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).tgt_spec{1}(1).value = 'mat';
            matlabbatch{14}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).tgt_spec{1}(2).name = 'strtype';
            matlabbatch{14}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).tgt_spec{1}(2).value = 'e';
            matlabbatch{14}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).sname = 'Liom Contrast Calculations: NIRS.mat';
            matlabbatch{14}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).src_exbranch = substruct('.','val', '{}',{13}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
            matlabbatch{14}.spm.tools.nirs10.model_estimate.liom_group.NIRSmat(1).src_output = substruct('.','NIRSmat');
            matlabbatch{14}.spm.tools.nirs10.model_estimate.liom_group.force_redo = 0;
            matlabbatch{14}.spm.tools.nirs10.model_estimate.liom_group.NIRSmatCopyChoice.NIRSmatOverwrite = struct([]);
            matlabbatch{14}.spm.tools.nirs10.model_estimate.liom_group.group_dir_name = 'GroupISA';
            matlabbatch{14}.spm.tools.nirs10.model_estimate.liom_group.number_dir_to_remove = 3;
            matlabbatch{14}.spm.tools.nirs10.model_estimate.liom_group.FFX_or_RFX = 1;
            matlabbatch{14}.spm.tools.nirs10.model_estimate.liom_group.ContrastChoice.automated_contrasts = struct([]);
            matlabbatch{14}.spm.tools.nirs10.model_estimate.liom_group.StatMethod = 1;
            matlabbatch{14}.spm.tools.nirs10.model_estimate.liom_group.contrast_p_value = 0.05;
            matlabbatch{14}.spm.tools.nirs10.model_estimate.liom_group.group_session_to_average = 1;
            matlabbatch{14}.spm.tools.nirs10.model_estimate.liom_group.simple_sum = 0;
            matlabbatch{14}.spm.tools.nirs10.model_estimate.liom_group.display_options.GenerateInverted = 1;
            matlabbatch{14}.spm.tools.nirs10.model_estimate.liom_group.display_options.GroupColorbars = 0;
            matlabbatch{14}.spm.tools.nirs10.model_estimate.liom_group.display_options.contrast_figures = 3;
            matlabbatch{14}.spm.tools.nirs10.model_estimate.liom_group.display_options.override_colorbar.colorbar_default = struct([]);
            matlabbatch{14}.spm.tools.nirs10.model_estimate.liom_group.display_options.figures_visible = 0;
            matlabbatch{14}.spm.tools.nirs10.model_estimate.liom_group.display_options.GroupFiguresIntoSubplots = 1;
            matlabbatch{14}.spm.tools.nirs10.model_estimate.liom_group.display_options.output_unc = 0;
            matlabbatch{14}.spm.tools.nirs10.model_estimate.liom_group.display_options.SmallFigures = 1;
            matlabbatch{14}.spm.tools.nirs10.model_estimate.liom_group.display_options.write_neg_pos = 0;
            matlabbatch{14}.spm.tools.nirs10.model_estimate.liom_group.display_options.save_nifti_contrasts = 0;
        end
        
    end
    save(fullfile(pathbatch,'Last_batch_main_module'),'matlabbatch');
    spm_jobman('run',matlabbatch);
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        disp(['Problem with ' cs]);
    end
end
