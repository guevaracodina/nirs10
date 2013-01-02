function coregnew1 = nirs_run_coreg_new_cfg
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'coreg');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration for coregistration: coreg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ForceReprocess      = cfg_menu;
ForceReprocess.tag  = 'ForceReprocess';
ForceReprocess.name = 'Force reprocessing of normalization';
ForceReprocess.labels = {'True','False'};
ForceReprocess.values = {1,0};
ForceReprocess.val  = {0};
ForceReprocess.help = {'Force reprocessing of the calculation of the normalization -- only.'}';

segT1_4fit         = cfg_files;
segT1_4fit.name    = 'Images segmented to fit P on.'; % The displayed name
segT1_4fit.tag     = 'segT1_4fit';
segT1_4fit.filter  = 'image';
segT1_4fit.num     = [0 Inf];
segT1_4fit.val{1}  = {''};
segT1_4fit.help    = {['(Optional) Choose the segmented image ',...
    '(0021_ or any other combination) to fit positions on scalp.',...
    'For multi-subject, the order of the images must correspond to ',...
    'the order of the subjects in the NIRS.mat matrix. ',...
    'If no input image is specified, the segmented image available in the ',...
    'NIRS matrix will be used (last one generated in MCsegment module).']};

fiducial_MNI_choice = cfg_menu;
fiducial_MNI_choice.tag    = 'fiducial_MNI_choice';
fiducial_MNI_choice.name   = 'Coordinate system for fiducials';
fiducial_MNI_choice.labels = {'Subject MNI','Normalized (i.e. template) MNI'};
fiducial_MNI_choice.values = {1,0};
fiducial_MNI_choice.val{1} = 0;
fiducial_MNI_choice.help   = {'Specify if coordinates of fiducials below are specified'
    'in the subject MNI coordinates or in the normalized MNI coordinate system. By default, the values'
    'are the fiducial positions in normalized MNI coordinates for the SPM8 standard subject.'}';

%Atlas Normalized coordinates of fiducials
nasion_wMNI         = cfg_entry; %nasion_wMNI
nasion_wMNI.name    = 'Nasion';
nasion_wMNI.tag     = 'nasion_wMNI';
nasion_wMNI.strtype = 'r';
nasion_wMNI.num     = [1 3];
nasion_wMNI.val{1}  = [0 84 -48];
nasion_wMNI.help    = {'Coordinates of the nasion.'};

AL_wMNI         = cfg_entry; %AL_wMNI
AL_wMNI.name    = 'Auricular left';
AL_wMNI.tag     = 'AL_wMNI';
AL_wMNI.strtype = 'r';
AL_wMNI.num     = [1 3];
AL_wMNI.val{1}  = [-83 -19 -38];
AL_wMNI.help    = {'Coordinates of Auricular Left.'};

AR_wMNI         = cfg_entry; %AR_wMNI
AR_wMNI.name    = 'Auricular right';
AR_wMNI.tag     = 'AR_wMNI';
AR_wMNI.strtype = 'r';
AR_wMNI.num     = [1 3];
AR_wMNI.val{1}  = [ 83 -19 -38];
AR_wMNI.help    = {'Coordinates of Auricular Right.'};

cortex_projection_method = nirs_dfg_cortex_projection_method;

render_template         = cfg_branch;
render_template.tag     = 'render_template';
render_template.name    = 'Render to SPM single subject template';
render_template.val     = {};
render_template.help    = {'Render to template.'};

render_subject         = cfg_branch;
render_subject.tag     = 'render_subject';
render_subject.name    = 'Render to subject';
render_subject.val     = {};
render_subject.help    = {'Render to subject.'}';

render_choice        = cfg_choice;
render_choice.name   = 'Projection (render) choice';
render_choice.tag    = 'render_choice';
render_choice.values = {render_template,render_subject};
render_choice.val    = {render_template};
render_choice.help   = {'Projection (render) choice: projection to template (normalized) or to subject'}';

coregnew1      = cfg_exbranch;
coregnew1.name = 'NIRScoreg (new)';
coregnew1.tag  = 'coregnew1';
coregnew1.val  = {NIRSmat redo1 ForceReprocess NIRSmatCopyChoice ...
    fiducial_MNI_choice nasion_wMNI AL_wMNI AR_wMNI ...
    cortex_projection_method render_choice};
coregnew1.prog = @nirs_run_coreg_new;
coregnew1.vout = @nirs_cfg_vout_coreg_new;
coregnew1.help = {'Automatic coregistration (new version).'};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_coreg_new(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
