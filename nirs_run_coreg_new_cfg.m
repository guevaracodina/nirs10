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

coreg_choice = cfg_menu;
coreg_choice.tag    = 'coreg_choice';
coreg_choice.name   = 'Coregistration Choice';
coreg_choice.labels = {'MNI Talairach Tournoux atlas','Subject Coordinates'};
coreg_choice.values = {1,0};
coreg_choice.val = {1};
coreg_choice.help   = {'Normalization choice. Need to specify appropriate'
    'file: normalized file if normalized to TT, unnormalized if in subject coordinates.'
    'Method in subject coordinates not coded up yet.'}';

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


render_file         = cfg_files;
render_file.name    = 'Render file';
render_file.tag     = 'render_file';
render_file.num     = [1 1];
render_file.help    = {'Grey matter (c1) anatomical image or rendered version of this c1 image.'
    'Normalized or not according to Normalization choice option below.'}';


render_normalize_choice = cfg_menu;
render_normalize_choice.tag    = 'render_normalize_choice';
render_normalize_choice.name   = 'Normalization Choice';
render_normalize_choice.labels = {'MNI Talairach Tournoux atlas','Subject Coordinates'};
render_normalize_choice.values = {1,0};
render_normalize_choice.val = {1};
render_normalize_choice.help   = {'Normalization choice. Need to specify appropriate'
    'file: normalized file if normalized to TT, unnormalized if in subject coordinates.'
    'Method in subject coordinates not coded up yet.'}';

render_template         = cfg_branch;
render_template.tag     = 'render_template';
render_template.name    = 'Render to SPM single subject template';
render_template.val     = {};
render_template.help    = {'Render to template.'};

occipital_shift         = cfg_entry;
occipital_shift.name    = 'Occipital shift';
occipital_shift.tag     = 'occipital_shift';
occipital_shift.strtype = 'r';
occipital_shift.num     = [1 1];
occipital_shift.val     = {0};
occipital_shift.help    = {'Enter shift in millimeters (positive to move more'
    'optodes toward the occipital view instead of on the frontal view)'
    'For a first try, 20 is suggested.'}';

render_subject         = cfg_branch;
render_subject.tag     = 'render_subject';
render_subject.name    = 'Render to subject';
render_subject.val     = {render_file render_normalize_choice occipital_shift};
render_subject.help    = {'Render to subject.'
    'Problem with coordinate systems and projections.'}';

render_choice        = cfg_choice;
render_choice.name   = 'Projection (render) choice';
render_choice.tag    = 'render_choice';
render_choice.values = {render_template,render_subject};
render_choice.val    = {render_template};
render_choice.help   = {'Projection (render) choice: projection to template (normalized), to subject'};

coregnew1      = cfg_exbranch;
coregnew1.name = 'NIRScoreg (new)';
coregnew1.tag  = 'coregnew1';
coregnew1.val  = {NIRSmat redo1 ForceReprocess NIRSmatCopyChoice ...
    fiducial_MNI_choice nasion_wMNI AL_wMNI AR_wMNI ...
    render_choice };
coregnew1.prog = @nirs_run_coreg;
coregnew1.vout = @nirs_cfg_vout_coreg;
coregnew1.help = {'Automatic coregistration (new version).'};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_coreg(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
