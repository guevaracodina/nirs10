function view3d1 = nirs_run_view3d_cfg
NIRSmat         = cfg_files; %Select NIRS.mat for this subject
NIRSmat.name    = 'NIRS.mat'; % The displayed name
NIRSmat.tag     = 'NIRSmat';       %file names
NIRSmat.filter  = 'mat';
NIRSmat.ufilter = '^NIRS.mat$';
NIRSmat.num     = [1 Inf];     % Number of inputs required
NIRSmat.help    = {'Select NIRS.mat for the subject(s).'}; % help text displayed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration for coregistration: view 3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

on_cortex        = cfg_menu;
on_cortex.tag    = 'on_cortex';
on_cortex.name   = 'View optodes on cortex.';
on_cortex.labels = {'Yes','No'};
on_cortex.values = {1,0};
on_cortex.def    = @(val)nirs_get_defaults('coregNIRS.view3d1.on_cortex', val{:});
on_cortex.help   = {'If answer is No, default view consists on displaying optodes on the scalp.'};

save_figure        = cfg_menu;
save_figure.tag    = 'save_figure';
save_figure.name   = 'Save 3D view figure.';
save_figure.labels = {'Yes','No'};
save_figure.values = {1,0};
save_figure.def    = @(val)nirs_get_defaults('coregNIRS.view3d1.save_figure', val{:});
save_figure.help   = {'Yes: 3D view will be saved as Matlab figure; No: figure will remain open.'};

segT1_4fit         = cfg_files;
segT1_4fit.name    = 'Images segmented viewed in 3D.'; % The displayed name
segT1_4fit.tag     = 'segT1_4fit';
segT1_4fit.filter  = 'image';
segT1_4fit.num     = [0 Inf];
segT1_4fit.val{1}  = {''};
segT1_4fit.help    = {['Optodes will be represented on the cortex or the ',...
    'scalp depending on the previous answer. If you answered Yes, please ',...
    'choose ''c3_'' image, otherwise choose the segmented image ',...
    '(0021_ or any other combination).',...
    'For multi-subject, the order of the images must correspond to ',...
    'the order of the subjects in the NIRS.mat matrix. ',...
    'If no input image is specified, the segmented image available in the ',...
    'NIRS matrix will be used (last one generated in MCsegment module).']};

% Executable Branch
view3d1      = cfg_exbranch;       % This is the branch that has information about how to run this module
view3d1.name = '3D view to check coregistration';             % The display name
view3d1.tag  = 'view3d1'; %Very important: tag is used when calling for execution
view3d1.val  = {NIRSmat on_cortex segT1_4fit save_figure};
view3d1.prog = @nirs_run_view3d;  % A function handle that will be called with the harvested job to run the computation
%view3d1.vout = @nirs_cfg_vout_view3d; % A function handle that will be called with the harvested job to determine virtual outputs
view3d1.help = {'Display 3D view of the skin with optodes so that user is able to check coregistration.'};

