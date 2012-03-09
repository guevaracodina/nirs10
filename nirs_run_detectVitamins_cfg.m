function detectVitamins1 = nirs_run_detectVitamins_cfg
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'Vit');

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
detectVitamins1.val  = {NIRSmat redo1 NIRSmatCopyChoice anatT1 output_prefix_woVit};
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