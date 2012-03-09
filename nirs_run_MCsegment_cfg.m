function MCsegment1 = nirs_run_MCsegment_cfg
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(0,'Seg');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration for MC segmentation: MCsegment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

image_in         = cfg_files;
image_in.tag     = 'image_in';
image_in.name    = 'Anatomical Images (Optional)';
image_in.help    = {'Select images to be processed. Same ROI will be selected ',...
    'for these images. Ancillary images will be saved... ',...
    'This module allows multi-subject processing, ',...
    'generating a segmented image for each subject. ',...
    'Note that a list of links to the segmented images will be ',...
    'available as a virtual output for further processing. ',...
    'It is essential that the order of the subjects be ',...
    'the same as in the readNIRS module.'};
image_in.filter  = 'image';
image_in.ufilter = '.*';
image_in.val{1}  = {''};
image_in.num     = [0 Inf];

%_______________________________________________________________________
output_autonaming      = cfg_menu;
output_autonaming.tag  = 'output_autonaming';
output_autonaming.name = 'Automatic output naming';
output_autonaming.labels = {'Yes','No'};
output_autonaming.values = {0,1};
output_autonaming.def = @(val)nirs_get_defaults('preprocANAT.MCsegment1.output_autonaming', val{:});
output_autonaming.help = {'Choose whether you want to choose the name of ',...
    'the output or not. If answer is ''Yes'', please enter name.'};
%_______________________________________________________________________
output_prefix      = cfg_entry;
output_prefix.tag  = 'output_prefix';
output_prefix.name = 'Prefix of the output image (MCsegment1)';
output_prefix.strtype = 's';
output_prefix.num     = [1 Inf];
output_prefix.def     = @(val)nirs_get_defaults('preprocANAT.MCsegment1.output_prefix', val{:});
output_prefix.help = {'You can choose to give a particular prefix to the ',...
    'output image. This prefix will be added at the left of the name of ',...
    'the image. A default name will be given otherwise.'};

%_______________________________________________________________________
sorting_method        = cfg_menu;
sorting_method.tag    = 'sorting_method';
sorting_method.name   = 'Sorting method';
sorting_method.labels = {'Median Filter','Opening','Gaussian Filter and Dilatation','Otsu','Median Filter and orthogonal Otsu'};
sorting_method.values = {0,1,2,3,4};
sorting_method.help   = {
    'NewSegment output images must be processed before pursuing segmentation : '
    '- (basic processing) Median Filter'
    '- Opening'
    '- Gaussian Filter and Dilatation'
    '- Otsu'
    'Explication de chacun des choix'
    }';

%_______________________________________________________________________
wtm      = cfg_branch;
wtm.tag  = 'wtm';
wtm.name = 'White matter';
sorting_method.def    = @(val)nirs_get_defaults('preprocANAT.MCsegment1.wtm', val{:});
wtm.val  = {sorting_method};
wtm.help = {'Options to process white matter image.'};
%_______________________________________________________________________
grm      = cfg_branch;
grm.tag  = 'grm';
grm.name = 'Grey matter';
sorting_method.def    = @(val)nirs_get_defaults('preprocANAT.MCsegment1.grm', val{:});
grm.val  = {sorting_method};
grm.help = {'Options to process grey matter image.'};
%_______________________________________________________________________
csf      = cfg_branch;
csf.tag  = 'csf';
csf.name = 'Cerebro-Spinal Fluid';
sorting_method.def    = @(val)nirs_get_defaults('preprocANAT.MCsegment1.csf', val{:});
csf.val  = {sorting_method};
csf.help = {'Options to process CSF image.'};
%_______________________________________________________________________
skl      = cfg_branch;
skl.tag  = 'skl';
skl.name = 'Skull';
sorting_method.def    = @(val)nirs_get_defaults('preprocANAT.MCsegment1.skl', val{:});
skl.val  = {sorting_method};
skl.help = {'Options to process skull image.'};
%_______________________________________________________________________
skn      = cfg_branch;
skn.tag  = 'skn';
skn.name = 'Skin';
sorting_method.def    = @(val)nirs_get_defaults('preprocANAT.MCsegment1.skn', val{:});
skn.val  = {sorting_method};
skn.help = {'Options to process skin image.'};

%_______________________________________________________________________
thresh_hs      = cfg_entry;
thresh_hs.tag  = 'thresh_hs';
thresh_hs.name = 'Threshold used to binarise head shadow in head shadow building';
thresh_hs.help = {'WHEN TO USE IT : if voxels are still zero after processing, you may increase this value keeping in mind, all boundaries will be blurred.'};
%_______________________________________________________________________
se_size_hs      = cfg_entry;
se_size_hs.tag  = 'se_size_hs';
se_size_hs.name = 'Size of the structural element in head shadow building';
se_size_hs.help = {
    'Size of the structural element (Mathematical morphology) used to clean the head shadow'
    'WHEN TO USE IT : (prefer change the threshold) nearly the same as the size of the threshold value.'
    }';
%_______________________________________________________________________
head_shadow      = cfg_branch;
head_shadow.tag  = 'head_shadow';
head_shadow.name = 'Build head shadowfor Monte Carlo Segmentation';

thresh_hs.def    = @(val)nirs_get_defaults('preprocANAT.head_shadow.thresh_hs', val{:});
se_size_hs.def   = @(val)nirs_get_defaults('preprocANAT.head_shadow.se_size_hs', val{:});
head_shadow.val  = {thresh_hs, se_size_hs};

head_shadow.help = {'The head shadow is obtained by summing all the ci-images (outputs from New Segment). It is used when sorting voxels inbetween one of the five layers. Once all the voxels have been sorted, some of them are zero whereas there are situated in the head volume, to solve this problem the head shadow is calculated and all the voxels belonging to the head are processed so that they are finally affected to one of the five layers'};

%_______________________________________________________________________
se_size_pi      = cfg_entry;
se_size_pi.tag  = 'se_size_pi';
se_size_pi.name = 'Size of the structural element in processing images';
se_size_pi.help = {
    'Size of the structural element (Mathematical morphology) used to clean images processed (ci-images ; outputs from New Segment)'
    'WHEN TO USE IT : .'
    }';
%_______________________________________________________________________
gaussfilt_size      = cfg_entry;
gaussfilt_size.tag  = 'gaussfilt_size';
gaussfilt_size.name = 'Size of the gaussian filter in processing images';
gaussfilt_size.help = {'...'
    'WHEN TO USE IT: if segmented image shows ''spatial instability'''
    '(regions with a lot of different layers)'}';
%_______________________________________________________________________
gaussfilt_sdev      = cfg_entry;
gaussfilt_sdev.tag  = 'gaussfilt_sdev';
gaussfilt_sdev.name = 'Standard deviation of the gaussian filter in processing images';
gaussfilt_sdev.help = {
    '...'
    'WHEN TO USE IT : '};
%_______________________________________________________________________
process_image      = cfg_branch;
process_image.tag  = 'process_image';
process_image.name = 'Process image for Monte Carlo Segmentation';

se_size_pi.def     = @(val)nirs_get_defaults('preprocANAT.process_image.se_size_pi', val{:});
gaussfilt_size.def = @(val)nirs_get_defaults('preprocANAT.process_image.gaussfilt_size', val{:});
gaussfilt_sdev.def = @(val)nirs_get_defaults('preprocANAT.process_image.gaussfilt_sdev', val{:});
process_image.val  = {se_size_pi, gaussfilt_size, gaussfilt_sdev};

process_image.help = {'Gaussian Filter is only used in method Gaussian Filter and dilate.'};

%_______________________________________________________________________
thresh_as      = cfg_entry;
thresh_as.tag  = 'thresh_as';
thresh_as.name = 'Threshold for voxels with no conflict of belonging';
thresh_as.def         = @(val)nirs_get_defaults('preprocANAT.MCsegment1.thresh_as', val{:});
thresh_as.help = {'Once all the ci-images have been calculated by SPM, ',...
    'the algorithm sorts voxels with respect to the number of layers ',...
    'they belong to. For those belonging to only one layer, the ',...
    'algorithm verifies the random variable is near enough from 1 ',...
    '(bigger than the said threshold).'};

%_______________________________________________________________________
vbm_seg      = cfg_menu;
vbm_seg.tag  = 'vbm_seg';
vbm_seg.name = 'Use VBM CSF and GM segmentation';
vbm_seg.labels = {'True','False'};
vbm_seg.values = {1,0};
vbm_seg.val  = {0};
vbm_seg.help = {'.'}';

%_______________________________________________________________________
rebel_surrounding      = cfg_entry;
rebel_surrounding.tag  = 'rebel_surrounding';
rebel_surrounding.name = 'Surrounding size for rebel voxels';
rebel_surrounding.def = @(val)nirs_get_defaults('preprocANAT.MCsegment1.rebel_surrounding', val{:});
rebel_surrounding.help = {'Size of the surrounding used to processed and hence sort rebel voxels.'};

%_______________________________________________________________________
rebel_thresh_hs      = cfg_entry;
rebel_thresh_hs.tag  = 'rebel_thresh_hs';
rebel_thresh_hs.name = 'Threshold for head shadow in Monte Carlo segmentation';
rebel_thresh_hs.def   = @(val)nirs_get_defaults('preprocANAT.MCsegment1.rebel_thresh_hs', val{:});
rebel_thresh_hs.help = {'Threshold to discriminate voxels in head shadow image.'};

force_reprocess      = cfg_menu;
force_reprocess.tag  = 'force_reprocess';
force_reprocess.name = 'Force reprocessing of segmentation';
force_reprocess.labels = {'True','False'};
force_reprocess.values = {1,0};
force_reprocess.val  = {0};
force_reprocess.help = {'Force reprocessing of segmentation.'
    'Note that SPM segmentation will never be reprocessed.'
    'This is only for reprocessing the preparation for Monte Carlo.'
    'Note also that previous results will be overwritten.'}';

MCsegment1      = cfg_exbranch;
MCsegment1.tag  = 'MCsegment1';
MCsegment1.name = 'MC Segmentation';
MCsegment1.val  = {NIRSmat redo1 NIRSmatCopyChoice force_reprocess  image_in output_autonaming ...
    output_prefix skn skl csf grm wtm vbm_seg thresh_as head_shadow ...
    rebel_surrounding rebel_thresh_hs process_image}; 
MCsegment1.prog = @nirs_run_MCsegment3;
MCsegment1.vout = @nirs_cfg_vout_MCsegment;
MCsegment1.help = {'Segmentation for Monte Carlo simulation and sorting ',...
    'of voxels of New Segment images. This module calls SPM''s New Segment ',...
    'unless it finds a c1 segmented file in the same directory as the input image.'};

%------------------------------------------------------------------------
function vout = nirs_cfg_vout_MCsegment(job)
% Determine what outputs will be present if this job is run. In this case,
% the structure of the inputs is fixed, and the output is always a single
% number. Note that input items may not be numbers, they can also be
% dependencies.
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
