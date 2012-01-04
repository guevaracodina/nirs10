function detectVitamins1 = nirs_run_detectVitamins_cfg

NIRSmat         = cfg_files; %Select NIRS.mat for this subject
NIRSmat.name    = 'NIRS.mat'; % The displayed name
NIRSmat.tag     = 'NIRSmat';       %file names
NIRSmat.filter  = 'mat';
NIRSmat.ufilter = '^NIRS.mat$';
NIRSmat.num     = [1 Inf];     % Number of inputs required
NIRSmat.help    = {'Select NIRS.mat for the subject(s).'}; % help text displayed

DelPreviousData      = cfg_menu;
DelPreviousData.tag  = 'DelPreviousData';
DelPreviousData.name = 'Delete Previous data file';
DelPreviousData.labels = {'True','False'};
DelPreviousData.values = {1,0};
DelPreviousData.val  = {0};
DelPreviousData.help = {'Delete the previous data file.'}';

CreateNIRSCopy_false         = cfg_branch;
CreateNIRSCopy_false.tag     = 'CreateNIRSCopy_false';
CreateNIRSCopy_false.name    = 'Do not copy NIRS structure';
CreateNIRSCopy_false.help    = {'Do not copy NIRS structure.'
    'This will write over the previous NIRS.mat'}';

NewNIRSdir         = cfg_entry;
NewNIRSdir.name    = 'Directory for NIRS.mat';
NewNIRSdir.tag     = 'NewNIRSdir';
NewNIRSdir.strtype = 's';
NewNIRSdir.val{1}    = 'NewDir';
NewNIRSdir.num     = [1 Inf];
NewNIRSdir.help    = {'Directory for NIRS.mat.'}';

CreateNIRSCopy         = cfg_branch;
CreateNIRSCopy.tag     = 'CreateNIRSCopy';
CreateNIRSCopy.name    = 'Create new directory and copy NIRS structure';
CreateNIRSCopy.val     = {NewNIRSdir};
CreateNIRSCopy.help    = {'Create new directory and copy NIRS structure there.'}';

%Common to most modules: for creating a new directory and copying NIRS.mat
NewDirCopyNIRS           = cfg_choice;
NewDirCopyNIRS.name      = 'Create new directory and copy NIRS.mat';
NewDirCopyNIRS.tag       = 'NewDirCopyNIRS';
NewDirCopyNIRS.values    = {CreateNIRSCopy_false CreateNIRSCopy};
NewDirCopyNIRS.val       = {CreateNIRSCopy_false};
NewDirCopyNIRS.help      = {'Choose whether to overwrite the NIRS.mat structure'
    'or to create a new directory'
    'and copy the NIRS.mat structure there'}';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration for preprocAnat: detect vitamins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

image_in         = cfg_files;
image_in.tag     = 'image_in';
image_in.name    = 'Image';
image_in.help    = {'Select the image to be processed. ROI will be selected in this image. Other stuff made with the image saved...'};
image_in.filter  = 'image';
image_in.ufilter = '.*';
image_in.num     = [1 1];

output_prefix_woVit         = cfg_entry;
output_prefix_woVit.name    = 'Prefix of the output anatomical image';
output_prefix_woVit.tag     = 'output_prefix_woVit';
output_prefix_woVit.strtype = 's';
output_prefix_woVit.num     = [1 Inf];
output_prefix_woVit.def  = @(val)nirs_get_defaults('preprocANAT.detectVitamins1.output_prefix_woVits', val{:});
output_prefix_woVit.help    = {'You can choose to give a particular prefix to the ',...
    'output image. This prefix will be added at the left of the name of the ',...
    'image. A default name will be given otherwise.'};

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

% Executable Branch
detectVitamins1      = cfg_exbranch;
detectVitamins1.name = 'Coregistration with fiducials';
detectVitamins1.tag  = 'detectVitamins1';
detectVitamins1.val  = {NIRSmat NewDirCopyNIRS anatT1 output_prefix_woVit};
detectVitamins1.prog = @nirs_run_detectVitamins;
detectVitamins1.vout = @nirs_cfg_vout_detectVitamins;
detectVitamins1.help = {['This module detects fiducial markers (vitamin capsules)'...
    ' on anatomical image (T1), saves their positions in the NIRS.mat matrix, and'...
    ' creates a copy of the anatomical image where the markers are erased.']};

% Make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_detectVitamins(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});