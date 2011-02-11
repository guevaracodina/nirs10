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
h13  = uimenu(h1,...
	'Label',	'Read LOT data',...
	'Separator',	'off',...
	'Tag',		'Read LOT data',...
	'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.readNIRS.lot1'');',...
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
h02  = uimenu(h00,...
    'Label',	'Permute Onsets',...
    'Separator',	'on',...
    'Tag',		'Permute Onsets',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.readOnsets.permuteOnsets'');',...
    'HandleVisibility','on');
h03  = uimenu(h00,...
    'Label',	'Add stimuli with HRFs for testing',...
    'Separator',	'on',...
    'Tag',		'Add stimuli with HRFs for testing',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.readOnsets.addTestStimuli'');',...
    'HandleVisibility','on');
h2  = uimenu(h0,...
    'Label',	 'Preprocess anatomical image',...
    'Separator', 'off',...
    'Tag',		 'Preprocess anatomical image',...
    'HandleVisibility','on');
h21  = uimenu(h2,...
    'Label',	'Segmentation for Monte Carlo',...
    'Separator',	'off',...
    'Tag',		'MCsegment',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.preprocANAT.MCsegment1'');',...
    'HandleVisibility','on');
% h22  = uimenu(h2,...
%     'Label',	'Remove vitamins from image',...
%     'Separator',	'off',...
%     'Tag',		'Remove vitamins from image',...
%     'HandleVisibility','on');
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
h32  = uimenu(h3,...
	'Label',	'Manual coregistration',...
	'Separator',	'off',...
	'Tag',		'Manual coregistration',...
	'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.coregNIRS.coreg_manual1'');',...
	'HandleVisibility','on');
% h33  = uimenu(h3,...
% 	'Label',	'Text coregistration',...
% 	'Separator',	'off',...
% 	'Tag',		'Text coregistration',...
% 	'HandleVisibility','on');
h34  = uimenu(h3,...
    'Label',	'Check (3D View)',...
    'Separator',	'on',...
    'Tag',		'Check (3D View)',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.coregNIRS.view3d1'');',...
    'HandleVisibility','on');
h35  = uimenu(h3,...
	'Label',	'Resize image',...
	'Separator',	'off',...
	'Tag',		'Resize image',...
	'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.coregNIRS.resize1'');',...
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
% h42 = uimenu(h4,...
% 	'Label',	'Preprocess: Mark Negative',...
% 	'Separator',	'off',...
% 	'Tag',		'Preprocess: Mark Negative',...
% 	'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.preprocessNIRS.mark_negative'');',...
% 	'HandleVisibility','on');
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
% h46  = uimenu(h4,...
% 	'Label',	'Filters: HPF and LPF',...
% 	'Separator',	'off',...
% 	'Tag',		'Filters: HPF and LPF',...
% 	'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.preprocessNIRS.HPF_LPF'');',...
% 	'HandleVisibility','on');
h47  = uimenu(h4,...
	'Label',	'Generate Analyzer files',...
	'Separator',	'on',...
	'Tag',		'Generate Analyzer files',...
	'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.preprocessNIRS.generate_vhdr_vmrk'');',...
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
h9  = uimenu(h0,...
    'Label',	'NIRS GL Model Specification',...
    'Separator',	'on',...
    'Tag',		'NIRS GL Model Specification',...
    'HandleVisibility','on');
h91  = uimenu(h9,...
    'Label',	'LIOM GLM Model Specification (WLS, BGLM, MDL (NIRS_SPM))',...
    'Separator',	'off',...
    'Tag',		'LIOM GLM Model Specification (WLS, BGLM, MDL (NIRS_SPM))',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.model_specify.wls_bglm_specify'');',...
    'HandleVisibility','on');
% h92  = uimenu(h9,...
%     'Label',	'NIRS_SPM (Korean) Model Specification',...
%     'Separator',	'off',...
%     'Tag',		'NIRS_SPM (Korean) Model Specification',...
%     'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.model_specify.NIRS_SPM_specify'');',...
%     'HandleVisibility','on');
% h93  = uimenu(h9,...
%     'Label',	'NIRS_SPM (Korean) Model Specification (NEW)',...
%     'Separator',	'off',...
%     'Tag',		'NIRS_SPM (Korean) Model Specification (NEW)',...
%     'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.model_specify.NIRS_SPM_specify_batch'');',...
%     'HandleVisibility','on');
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
% h100  = uimenu(h10,...
% 	'Label',	'NIRS_SPM Filters: HPF and LPF',...
% 	'Separator',	'off',...
% 	'Tag',		'NIRS_SPM Filters: HPF and LPF',...
% 	'CallBack','spm_jobman(''interactive'','''',''spm.tools.nirs10.model_estimate.NIRS_SPM_HPF_LPF'');',...
% 	'HandleVisibility','on');
% h102  = uimenu(h10,...
%     'Label',	'NIRS_SPM (Korean) Model Estimation',...
%     'Separator',	'off',...
%     'Tag',		'NIRS_SPM (Korean) Model Estimation',...
%     'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.model_estimate.NIRS_SPM_estimate'');',...
%     'HandleVisibility','on');
% h103  = uimenu(h10,...
%     'Label',	'NIRS_SPM (Korean) Model Estimation (NEW, batch)',...
%     'Separator',	'off',...
%     'Tag',		'NIRS_SPM (Korean) Model Estimation (NEW, batch)',...
%     'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.model_estimate.NIRS_SPM_estimate_batch'');',...
%     'HandleVisibility','on');
% h104  = uimenu(h10,...
%     'Label',	'NIRS_SPM (Korean) Contrast Calculation',...
%     'Separator',	'off',...
%     'Tag',		'NIRS_SPM (Korean) Contrast Calculation',...
%     'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.model_estimate.NIRS_SPM_contrast'');',...
%     'HandleVisibility','on');
h104b  = uimenu(h10,...
    'Label',	'Liom Contrast Calculation',...
    'Separator',	'off',...
    'Tag',		'Liom Contrast Calculation',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.model_estimate.liom_contrast'');',...
    'HandleVisibility','on');
% h105  = uimenu(h10,...
%     'Label',	'NIRS_SPM (Korean) Group Model Estimation',...
%     'Separator',	'off',...
%     'Tag',		'NIRS_SPM (Korean) Group Model Estimation',...
%     'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.model_estimate.NIRS_SPM_group'');',...
%     'HandleVisibility','on');
h105b  = uimenu(h10,...
    'Label',	'Liom Group Model Estimation',...
    'Separator',	'off',...
    'Tag',		'Liom Group Model Estimation',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.model_estimate.liom_group'');',...
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
h11  = uimenu(h0,...
    'Label',	'NIRS GLM Model Results Display',...
    'Separator',	'off',...
    'Tag',		'NIRS GLM Model Results Display',...
    'HandleVisibility','on');
h111  = uimenu(h11,...
    'Label',	'NIRS_SPM GLM Model Results Display',...
    'Separator',	'off',...
    'Tag',		'NIRS_SPM GLM Model Results Display',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.model_display.NIRS_SPM_model_display'');',...
    'HandleVisibility','on');
h112  = uimenu(h11,...
    'Label',	'NIRS_SPM Contrast Estimates Display',...
    'Separator',	'off',...
    'Tag',		'NIRS_SPM Contrast Estimates Display',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.model_display.NIRS_SPM_contrast_display'');',...
    'HandleVisibility','on');
h113  = uimenu(h11,...
    'Label',	'NIRS_SPM Detrending and Protocol Diagnostic',...
    'Separator',	'off',...
    'Tag',		'NIRS_SPM Detrending and Protocol Diagnostic',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.model_display.NIRS_SPM_diagnostic'');',...
    'HandleVisibility','on');
h12  = uimenu(h0,...
    'Label',	'NIRS HDM',...
    'Separator',	'off',...
    'Tag',		'NIRS HDM',...
    'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.NIRS_HDM'');',...
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
% h14  = uimenu(h0,...
%     'Label',	'NIRS DCM',...
%     'Separator',	'off',...
%     'Tag',		'NIRS DCM',...
%     'CallBack',  'spm_jobman(''interactive'','''',''spm.tools.nirs10.NIRS_DCM'');',...
%     'HandleVisibility','on');


