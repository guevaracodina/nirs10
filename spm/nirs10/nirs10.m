function nirs10
% Based on spm_vbm8 from vbm8
% Toolbox wrapper to call mctools functions

rev = '$Rev: 1 $';

SPMid = spm('FnBanner',mfilename,rev);
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','NIRS10');
spm_help('!ContextHelp',mfilename);
spm_help('!Disp','nirs10.man','',Fgraph,'NIRS10');

fig = spm_figure('GetWin','Interactive');
h0  = uimenu(fig,...
    'Label',	'nirs10',...
    'Separator',	'on',...
    'Tag',		'NIRS',...
    'HandleVisibility','on');
h1  = uimenu(h0,...
	'Label',	'Read NIRS data',...
	'Separator',	'on',...
	'Tag',		'Read NIRS data',...
	'HandleVisibility','on');
h11  = uimenu(h1,...
	'Label',	'Read BOXY data',...
	'Separator',	'off',...
	'Tag',		'Read BOXY data',...
	'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.readNIRS.boxy1'');',...
	'HandleVisibility','on');
h12  = uimenu(h1,...
	'Label',	'Read and format CRIUGM data',...
	'Separator',	'off',...
	'Tag',		'Read and format CRIUGM data',...
	'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.readNIRS.criugm1'');',...
	'HandleVisibility','on');
h11b  = uimenu(h1,...
	'Label',	'Read BOXY data - Manual',...
	'Separator',	'off',...
	'Tag',		'Read BOXY data - Manual',...
	'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.readNIRS.boxy_manual1'');',...
	'HandleVisibility','on');
h00  = uimenu(h0,...
	'Label',	'Read Onsets',...
	'Separator',	'off',...
	'Tag',		'Read Onsets',...
	'HandleVisibility','on');
h01  = uimenu(h00,...
    'Label',	'Read NIRS Analyzer 2 Onsets',...
    'Separator',	'off',...
    'Tag',		'Read NIRS Analyzer 2 Onsets',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.readOnsets.AnalyzerOnsets'');',...
    'HandleVisibility','on');
h01b  = uimenu(h00,...
    'Label',	'Read Eprime onsets',...
    'Separator',	'off',...
    'Tag',		'Read Eprime onsets',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.readOnsets.readEprimeOnsets1'');',...
    'HandleVisibility','on');
h02  = uimenu(h00,...
    'Label',	'Permute Onsets',...
    'Separator',	'on',...
    'Tag',		'Permute Onsets',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.readOnsets.permuteOnsets1'');',...
    'HandleVisibility','on');
h03  = uimenu(h00,...
    'Label',	'Add stimuli with HRFs for testing',...
    'Separator',	'on',...
    'Tag',		'Add stimuli with HRFs for testing',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.readOnsets.addTestStimuli'');',...
    'HandleVisibility','on');
h4  = uimenu(h0,...
	'Label',	'NIRS Preprocessing',...
	'Separator',	'off',...
	'Tag',		'NIRS Preprocessing',...
	'HandleVisibility','on');
h40 = uimenu(h4,...
	'Label',	'Preprocess: Remove Channels (stdev)',...
	'Separator',	'off',...
	'Tag',		'Preprocess: Remove Channels (stdev)',...
	'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.preprocessNIRS.remove_chn_stdev'');',...
	'HandleVisibility','on');
h41 = uimenu(h4,...
	'Label',	'Preprocess: Heart Rate Utility',...
	'Separator',	'off',...
	'Tag',		'Preprocess: Heart Rate Utility',...
	'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.preprocessNIRS.criugm_paces1'');',...
	'HandleVisibility','on');
h43 = uimenu(h4,...
	'Label',	'Preprocess: Mark Movement',...
	'Separator',	'off',...
	'Tag',		'Preprocess: Mark Movement',...
	'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.preprocessNIRS.mark_movement'');',...
	'HandleVisibility','on');
h44  = uimenu(h4,...
	'Label',	'Normalize to baseline',...
	'Separator',	'off',...
	'Tag',		'Normalize to baseline',...
	'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.preprocessNIRS.normalize_baseline'');',...
	'HandleVisibility','on');
h45  = uimenu(h4,...
	'Label',	'Convert OD to HbO/HbR',...
	'Separator',	'off',...
	'Tag',		'ODtoHbOHbR',...
	'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.preprocessNIRS.ODtoHbOHbR'');',...
	'HandleVisibility','on');
h47  = uimenu(h4,...
	'Label',	'Generate Analyzer files',...
	'Separator',	'on',...
	'Tag',		'Generate Analyzer files',...
	'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.preprocessNIRS.generate_vhdr_vmrk'');',...
	'HandleVisibility','on');
h2  = uimenu(h0,...
    'Label',	 'Preprocess anatomical image',...
    'Separator', 'off',...
    'Tag',		 'Preprocess anatomical image',...
    'HandleVisibility','on');
h20  = uimenu(h2,...
    'Label',	'Detect fiducials in anatomical image',...
    'Separator',	'off',...
    'Tag',		'Detect fiducials',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.preprocANAT.detectVitamins1'');',...
    'HandleVisibility','on');
h21  = uimenu(h2,...
    'Label',	'Segmentation for Monte Carlo',...
    'Separator',	'off',...
    'Tag',		'MCsegment',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.preprocANAT.MCsegment1'');',...
    'HandleVisibility','on');
h23  = uimenu(h2,...
    'Label',	'Select Region of Interest',...
    'Separator',	'off',...
    'Tag',		'Select ROI',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.preprocANAT.buildroi1'');',...
    'HandleVisibility','on');
h3  = uimenu(h0,...
	'Label',	'Coregister NIRS data',...
	'Separator',	'off',...
	'Tag',		'Coregister NIRS data',...
	'HandleVisibility','on');
h31  = uimenu(h3,...
	'Label',	'Automatic coregistration',...
	'Separator',	'off',...
	'Tag',		'Automatic coregistration',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.coregNIRS.coreg1'');',...	
	'HandleVisibility','on');
h31new  = uimenu(h3,...
	'Label',	'Automatic coregistration (new)',...
	'Separator',	'off',...
	'Tag',		'Automatic coregistration (new)',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.coregNIRS.coregnew1'');',...	
	'HandleVisibility','on');
h31bis  = uimenu(h3,...
	'Label',	'Automatic coregistration to template',...
	'Separator',	'off',...
	'Tag',		'Automatic coregistration to template',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.coregNIRS.coreg2'');',...	
	'HandleVisibility','on');
h31ter  = uimenu(h3,...
	'Label',	'Automatic coregistration of template helmet',...
	'Separator',	'off',...
	'Tag',		'Automatic coregistration of template helmet',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.coregNIRS.coreg3'');',...	
	'HandleVisibility','on');
h32  = uimenu(h3,...
	'Label',	'Manual coregistration',...
	'Separator',	'off',...
	'Tag',		'Manual coregistration',...
	'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.coregNIRS.coreg_manual1'');',...
	'HandleVisibility','on');
h34  = uimenu(h3,...
    'Label',	'Check (3D View)',...
    'Separator',	'on',...
    'Tag',		'Check (3D View)',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.coregNIRS.view3d1'');',...
    'HandleVisibility','on');
h5  = uimenu(h0,...
	'Label',	'Configure Monte Carlo inputs',...
	'Separator',	'on',...
	'Tag',		'Configure Monte Carlo inputs',...
	'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.configMC1'');',...
	'HandleVisibility','on');
h6  = uimenu(h0,...
    'Label',	'Run Monte Carlo Simulation',...
    'Separator',	'off',...
    'Tag',		'Run Monte Carlo Simulation',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.runMC1'');',...
    'HandleVisibility','on');
h7  = uimenu(h0,...
    'Label',	'Assemble sensitivity matrix',...
    'Separator',	'off',...
    'Tag',		'Assemble sensitivity matrix',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.makesens1'');',...
    'HandleVisibility','on');
h14  = uimenu(h0,...
	'Label',	'Calculate Partial Volume Effect',...
	'Separator',	'off',...
	'Tag',		'Calculate Partial Volume Effect',...
	'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.calculatePVE1'');',...
	'HandleVisibility','on');
h8  = uimenu(h0,...
    'Label',	'3D Reconstruction of NIRS data',...
    'Separator',	'off',...
    'Tag',		'3D Reconstruction of NIRS data',...
    'HandleVisibility','on');
h81  = uimenu(h8,...
    'Label',	'Tikhonov',...
    'Separator',	'off',...
    'Tag',		'Tikhonov',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.model_reconstruct.tikhonov1'');',...
    'HandleVisibility','on');
h82  = uimenu(h8,...
    'Label',	'ReML',...
    'Separator',	'off',...
    'Tag',		'ReML',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.model_reconstruct.ReMLreconstruct1'');',...
    'HandleVisibility','on');
h83  = uimenu(h8,...
    'Label',	'Test reconstruction with phantom',...
    'Separator',	'on',...
    'Tag',		'Test reconstruction with phantom',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.model_reconstruct.testreconstruct1'');',...
    'HandleVisibility','on');
h84  = uimenu(h8,...
    'Label',	'Check reconstruction',...
    'Separator',	'off',...
    'Tag',		'Check reconstruction',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.model_reconstruct.checkreconstruct1'');',...
    'HandleVisibility','on');
h9  = uimenu(h0,...
    'Label',	'NIRS GL Model Specification',...
    'Separator',	'on',...
    'Tag',		'NIRS GL Model Specification',...
    'HandleVisibility','on');
h90  = uimenu(h9,...
    'Label',	'LIOM Intrasubject Average',...
    'Separator',	'off',...
    'Tag',		'LIOM Intrasubject Average',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.model_specify.liom_intrasubject_average'');',...
    'HandleVisibility','on');
h91  = uimenu(h9,...
    'Label',	'LIOM GLM Model Specification (WLS, BGLM, MDL (NIRS_SPM))',...
    'Separator',	'off',...
    'Tag',		'LIOM GLM Model Specification (WLS, BGLM, MDL (NIRS_SPM))',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.model_specify.wls_bglm_specify'');',...
    'HandleVisibility','on');
h92  = uimenu(h9,...
    'Label',	'LIOM CINE',...
    'Separator',	'off',...
    'Tag',		'LIOM CINE',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.model_specify.liom_cine'');',...
    'HandleVisibility','on');
h10  = uimenu(h0,...
    'Label',	'NIRS GL Model Estimation',...
    'Separator',	'off',...
    'Tag',		'NIRS GL Model Estimation',...
    'HandleVisibility','on');
h101  = uimenu(h10,...
    'Label',	'Liom GLM Estimation (WLS, BGLM, MDL (NIRS_SPM))',...
    'Separator',	'off',...
    'Tag',		'Liom GLM Estimation (WLS, BGLM, MDL (NIRS_SPM))',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.model_estimate.wls_bglm_estimate'');',...
    'HandleVisibility','on');
h104b  = uimenu(h10,...
    'Label',	'Liom Contrast Calculation',...
    'Separator',	'off',...
    'Tag',		'Liom Contrast Calculation',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.model_estimate.liom_contrast'');',...
    'HandleVisibility','on');
h105b  = uimenu(h10,...
    'Label',	'Liom Group Model Estimation',...
    'Separator',	'off',...
    'Tag',		'Liom Group Model Estimation',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.model_estimate.liom_group'');',...
    'HandleVisibility','on');
h105c = uimenu(h10,...
    'Label',	'Extract map data',...
    'Separator',	'off',...
    'Tag',		'Extract map data',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.model_estimate.extract_map_data'');',...
    'HandleVisibility','on');
h105cbis = uimenu(h10,...
    'Label',	'Extract map data -- simplified',...
    'Separator',	'off',...
    'Tag',		'Extract map data -- simplified',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.model_estimate.extract_map_data_simplified'');',...
    'HandleVisibility','on');
h105d = uimenu(h10,...
    'Label',	'One-way anova',...
    'Separator',	'off',...
    'Tag',		'One-way anova',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.model_estimate.liom_1way_anova'');',...
    'HandleVisibility','on');
h105e = uimenu(h10,...
    'Label',	'Two-way anova',...
    'Separator',	'off',...
    'Tag',		'Two-way anova',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.model_estimate.liom_2way_anova'');',...
    'HandleVisibility','on');
h105f = uimenu(h10,...
    'Label',	'Mixed two-way anova',...
    'Separator',	'off',...
    'Tag',		'Mixed two-way anova',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.model_estimate.liom_mixed2way_anova'');',...
    'HandleVisibility','on');
h106a  = uimenu(h10,...
    'Label',	'Analyze GLMs',...
    'Separator',	'on',...
    'Tag',		'Analyze GLMs',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.model_estimate.AnalyzeGLM'');',...
    'HandleVisibility','on');
h106  = uimenu(h10,...
    'Label',	'ROC - sensitivity and specificity',...
    'Separator',	'on',...
    'Tag',		'ROC - sensitivity and specificity',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.model_estimate.ROCtest'');',...
    'HandleVisibility','on');
h12  = uimenu(h0,...
    'Label',	'Hemodynamic modelling',...
    'Separator',	'off',...
    'Tag',		'Hemodynamic modelling',...
    'HandleVisibility','on');
h121  = uimenu(h12,...
    'Label',	'HDM (new)',...
    'Separator',	'off',...
    'Tag',		'HDM (new)',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.HDM_options.HDM'');',...
    'HandleVisibility','on');
h122  = uimenu(h12,...
    'Label',	'NIRS HDM (old)',...
    'Separator',	'off',...
    'Tag',		'NIRS HDM (old)',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.HDM_options.NIRS_HDM'');',...
    'HandleVisibility','on');
h123  = uimenu(h12,...
    'Label',	'LIOM HDM (old, for BOLD)',...
    'Separator',	'off',...
    'Tag',		'LIOM HDM (old, for BOLD)',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.HDM_options.liom_HDM'');',...
    'HandleVisibility','on');
h124  = uimenu(h12,...
    'Label',	'SCKS',...
    'Separator',	'off',...
    'Tag',		'SCKS',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.HDM_options.SCKS'');',...
    'HandleVisibility','on');
h13 = uimenu(h0,...
    'Label',	'CRIUGM',...
    'Separator',	'on',...
    'Tag',		'CRIUGM',...
    'HandleVisibility','on');
h131 = uimenu(h13,...
    'Label',	'VOIRE project',...
    'Separator',	'off',...
    'Tag',		'VOIRE project',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.CRIUGM.runVOIRE1'');',...
    'HandleVisibility','on');
h131 = uimenu(h13,...
    'Label',	'Mobility project',...
    'Separator',	'off',...
    'Tag',		'Mobility project',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.CRIUGM.runMOB1'');',...
    'HandleVisibility','on');

% ------------------------------------------------------------------------------
% functional connectivity NIRS (fcNIRS) menus
% ------------------------------------------------------------------------------
h13a = uimenu(h0,...
    'Label',            'fcNIRS',...
    'Separator',        'on',...
    'Tag',              'fcNIRS_menu',...
    'HandleVisibility', 'on');
h13a1 = uimenu(h139,...
    'Label',            'BPF & Downsampling',...
    'Separator',        'off',...
    'Tag',              'BPF & Downsampling',...
    'CallBack',         'spm_jobman(''interactive'','''',''spm.tools.nirs10.fcNIRS.filtdown1'');',...
    'HandleVisibility', 'on');
h13a2 = uimenu(h139,...
    'Label',            'Send e-mail',...
    'Separator',        'off',...
    'Tag',              'Send e-mail',...
    'CallBack',         'spm_jobman(''interactive'','''',''spm.tools.nirs10.fcNIRS.sendmail1'');',...
    'HandleVisibility', 'on');
% ------------------------------------------------------------------------------

h16  = uimenu(h0,...
    'Label',	'NIRS Utilities',...
    'Separator',	'on',...
    'Tag',		'NIRS Utilities',...
    'HandleVisibility','on');
h161  = uimenu(h16,...
    'Label',	'Convert 2D topo nii to 2D topo views',...
    'Separator',	'off',...
    'Tag',		'Convert 2D topo nii to 2D topo views',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.nirs_utilities.convert_nii_to_2Dtopo'');',...
    'HandleVisibility','on');
h162  = uimenu(h16,...
	'Label',	'Projection of focus on cortex',...
	'Separator',	'off',...
	'Tag',		'Projection of focus on cortex',...
	'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.nirs_utilities.liom_projection'');',...
	'HandleVisibility','on');
h163  = uimenu(h16,...
	'Label',	'Sensitivity and specificity check',...
	'Separator',	'off',...
	'Tag',		'Sensitivity and specificity check',...
	'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.nirs_utilities.liom_sensitivity_specificity'');',...
	'HandleVisibility','on');
h164  = uimenu(h16,...
	'Label',	'Coregistration on orthogonal image',...
	'Separator',	'off',...
	'Tag',		'Coregistration on orthogonal image',...
	'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.nirs_utilities.liom_OrthCoreg'');',...
	'HandleVisibility','on');
h17  = uimenu(h0,...
    'Label',	'Anatomical MRI Utilities',...
    'Separator',	'on',...
    'Tag',		'Anatomical MRI Utilities',...
    'HandleVisibility','on');
h171  = uimenu(h17,...
    'Label',	'Thickness below selected points',...
    'Separator',	'off',...
    'Tag',		'Thickness below selected points',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.aMRI_utilities.Uthickness'');',...
    'HandleVisibility','on');
h172  = uimenu(h17,...
	'Label',	'Resize image',...
	'Separator',	'off',...
	'Tag',		'Resize image',...
	'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.aMRI_utilities.resize1'');',...
	'HandleVisibility','on');
h140 = uimenu(h0,...
    'Label',	'Full Criugm study',...
    'Separator',	'off',...
    'Tag',		'Full Criugm study',...
    'CallBack',  'spm_jobman(''interactive'',''fullStudyCriugm_topo_cfg.m'');',...
    'HandleVisibility','on');
% h150 = uimenu(h0,...
%     'Label',	'Full epilepsy study',...
%     'Separator',	'off',...
%     'Tag',		'Full epilepsy study',...
%     'CallBack',  'spm_jobman(''interactive'',''fullStudyEpilepsy_topo_cfg.m'');',...
%     'HandleVisibility','on');
h151 = uimenu(h0,...
    'Label',	'Full epilepsy study - Flexible Select',...
    'Separator',	'off',...
    'Tag',		'Full epilepsy study - Flexible Select',...
    'CallBack',  'spm_jobman(''interactive'',''fullStudyEpilepsyManual_topo_cfg.m'');',...
    'HandleVisibility','on');
h152 = uimenu(h0,...
    'Label',	'GLM + Contrasts + Group',...
    'Separator',	'off',...
    'Tag',		'GLM + Contrasts + Group',...
    'CallBack',  'spm_jobman(''interactive'',''GLM_contrasts_group_cfg.m'');',...
    'HandleVisibility','on');
