function coreg1 = nirs_run_coreg_cfg
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'coreg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration for coregistration: coreg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

template_mode = cfg_menu;
template_mode.tag    = 'template_mode';
template_mode.name   = 'Using a template for all subjects?';
template_mode.labels = {'Yes','No'};
template_mode.values = {1,0};
template_mode.val{1} = 0;
template_mode.help   = {'If coregistration is to be performed for each subject '
    'separately, choose No. If coregistration should be performed only once, '
    'for the first subject, since a tem47plate has been specified to use for all subjects, choose Yes.'}';

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

fid_in_subject_MNI = cfg_menu;
fid_in_subject_MNI.tag    = 'fid_in_subject_MNI';
fid_in_subject_MNI.name   = 'Fiducials in subject MNI coordinates?';
fid_in_subject_MNI.labels = {'Yes','No'};
fid_in_subject_MNI.values = {1,0};
fid_in_subject_MNI.def    = @(val)nirs_get_defaults('coregNIRS.coreg1.fid_in_subject_MNI', val{:});
fid_in_subject_MNI.help   = {'Specify if coordinates of fiducials below are specified'
    'in the subject or patient MNI coordinates. The default option is NO: the default values'
    'are fiducial positions in normalized MNI coordinates of the SPM8 standard subject.'}';

%Atlas Normalized coordinates of fiducials
nasion_wMNI         = cfg_entry; %nasion_wMNI
nasion_wMNI.name    = 'Nasion';
nasion_wMNI.tag     = 'nasion_wMNI';
nasion_wMNI.strtype = 'r';
nasion_wMNI.num     = [1 3];
nasion_wMNI.def     = @(val)nirs_get_defaults('coregNIRS.coreg1.nasion_wMNI', val{:});
nasion_wMNI.help    = {'Coordinates of the nasion.'};

AL_wMNI         = cfg_entry; %AL_wMNI
AL_wMNI.name    = 'Auricular left';
AL_wMNI.tag     = 'AL_wMNI';
AL_wMNI.strtype = 'r';
AL_wMNI.num     = [1 3];
AL_wMNI.def     = @(val)nirs_get_defaults('coregNIRS.coreg1.AL_wMNI', val{:});
AL_wMNI.help    = {'Coordinates of Auricular Left.'};

AR_wMNI         = cfg_entry; %AR_wMNI
AR_wMNI.name    = 'Auricular right';
AR_wMNI.tag     = 'AR_wMNI';
AR_wMNI.strtype = 'r';
AR_wMNI.num     = [1 3];
AR_wMNI.def     = @(val)nirs_get_defaults('coregNIRS.coreg1.AR_wMNI', val{:});
AR_wMNI.help    = {'Coordinates of Auricular Right.'};

GenDataTopo = cfg_menu;
GenDataTopo.tag    = 'GenDataTopo';
GenDataTopo.name   = 'Generate topo data.';
GenDataTopo.labels = {'Yes','No'};
GenDataTopo.values = {1,0};
GenDataTopo.def    = @(val)nirs_get_defaults('coregNIRS.coreg1.GenDataTopo', val{:});
GenDataTopo.help   = {'Generate rend data (NIRS_SPM) for topographic '
    'reconstruction - stored in a separate file: TopoData.mat'}';

View6Projections      = cfg_menu;
View6Projections.tag  = 'View6Projections';
View6Projections.name = 'View the 6 Projections';
View6Projections.labels = {'True','False'};
View6Projections.values = {1,0};
View6Projections.val  = {0};
View6Projections.help = {'View channel locations for the 6 projections.'}';

Save6Projections      = cfg_menu;
Save6Projections.tag  = 'Save6Projections';
Save6Projections.name = 'Save the 6 Projections';
Save6Projections.labels = {'True','False'};
Save6Projections.values = {1,0};
Save6Projections.val  = {1};
Save6Projections.help = {'Save images of channel locations for the 6 projections.'}';

ForceReprocess      = cfg_menu;
ForceReprocess.tag  = 'ForceReprocess';
ForceReprocess.name = 'Force reprocessing of normalization';
ForceReprocess.labels = {'True','False'};
ForceReprocess.values = {1,0};
ForceReprocess.val  = {0};
ForceReprocess.help = {'Force reprocessing of the calculation of the normalization -- only.'}';

render_file         = cfg_files;
render_file.name    = 'Render file';
render_file.tag     = 'render_file';
%render_file.filter  = 'image';
%render_file.ufilter = '.*';
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
render_subject.help    = {'Render to subject. OPTION NOT FUNCTIONAL YET: '
    'Problem with coordinate systems and projections.'}';

render_choice        = cfg_choice;
render_choice.name   = 'Render to template or subject';
render_choice.tag    = 'render_choice';
render_choice.values = {render_template,render_subject};
render_choice.val    = {render_template};
render_choice.help   = {'Render to template or subject.'};

anatT1         = cfg_files; %Select T1 for this subject
anatT1.name    = 'Raw anatomical image (optional)'; % The displayed name
anatT1.tag     = 'anatT1';       %file names
anatT1.filter  = 'image';
anatT1.ufilter = '.*';
anatT1.val{1}  = {''};
anatT1.num     = [0 Inf];     % Number of inputs required
anatT1.help    = {'CORRECTION: THIS IS NO LONGER OPTIONAL.'
    'It is required for the coregistration module and the contrast module'
    'Optional, can be specified in MC Segment, or earlier '
    'and be available in NIRS.mat structure.'
    'Select raw anatomical image(s) for the subject(s). '
    'If several subjects, the images '
    'must be in the same order as the NIRS.mat structures.'}';

anatT1_template         = cfg_files;
anatT1_template.name    = 'Anatomical template image';
anatT1_template.tag     = 'anatT1_template';       %file names
anatT1_template.filter  = 'image';
anatT1_template.ufilter = '.*';
anatT1_template.def = @(val)nirs_get_defaults('coregNIRS.coreg1.anatT1_template', val{:});
anatT1_template.num     = [1 1];     % Number of inputs required
anatT1_template.help    = {'Select anatomical template image for this subject.'};

coreg1      = cfg_exbranch;
coreg1.name = 'NIRScoreg';
coreg1.tag  = 'coreg1';
coreg1.val  = {NIRSmat redo1 NIRSmatCopyChoice template_mode anatT1 segT1_4fit ...
    anatT1_template fid_in_subject_MNI nasion_wMNI AL_wMNI AR_wMNI ...
    GenDataTopo render_choice View6Projections Save6Projections ForceReprocess};
coreg1.prog = @nirs_run_coreg;
coreg1.vout = @nirs_cfg_vout_coreg;
coreg1.help = {'Automatic coregistration.'};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_coreg(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
