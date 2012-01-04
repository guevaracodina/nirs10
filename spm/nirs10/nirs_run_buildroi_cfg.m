function buildroi1 = nirs_run_buildroi_cfg
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
%Configuration for MC segmentation: buildroi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

keepChannels         = cfg_entry;
keepChannels.name    = 'List of channels to keep';
keepChannels.tag     = 'keepChannels';
keepChannels.strtype = 'r';
keepChannels.num     = [1 Inf];
keepChannels.def     = @(val)nirs_get_defaults('readOnsets.addTestStimuli.keepChannels', val{:});
keepChannels.help    = {'Enter channel numbers to keep.'};

AllChannels           = cfg_branch;
AllChannels.name      = 'Keep ALL channels';
AllChannels.tag       = 'AllChannels';
AllChannels.help      = {'Keep all channels'};

keepAllChannels        = cfg_choice;
keepAllChannels.name   = 'Specify channels to keep';
keepAllChannels.tag    = 'keepAllChannels';
keepAllChannels.values = {AllChannels keepChannels};
%Do not know how to specify a default value by a call using .def for a
%cfg_choice object
keepAllChannels.val    = {keepChannels};
keepAllChannels.help   = {'Choose whether to keep all channels or select a subset.'};

image_in         = cfg_files;
image_in.tag     = 'image_in';
image_in.name    = 'Image';
image_in.help    = {'Select the image to be processed. ROI will be selected in this image. Other stuff made with the image saved...'};
image_in.filter  = 'image';
image_in.ufilter = '.*';
image_in.val{1} = {''};
image_in.num     = [0 1];

%_______________________________________________________________________
output_prefix         = cfg_entry;
output_prefix.name    = 'Prefix of the output image (buildroi1)';
output_prefix.tag     = 'output_prefix';
output_prefix.strtype = 's';
output_prefix.num     = [1 Inf];
output_prefix.def  = @(val)nirs_get_defaults('preprocANAT.buildroi1.output_prefix', val{:});
output_prefix.help    = {'You can choose to give a particular prefix to the ',...
    'output image. This prefix will be added at the left of the name of the ',...
    'image. A default name will be given otherwise.'};

%_______________________________________________________________________
buildroi1      = cfg_exbranch;
buildroi1.tag  = 'buildroi1';
buildroi1.name = 'Build ROI';
buildroi1.val  = {NIRSmat NewDirCopyNIRS keepAllChannels image_in output_prefix};
buildroi1.prog = @nirs_run_buildroi2;
buildroi1.vout = @nirs_cfg_vout_buildroi1;
buildroi1.help = {'Define region of interest containing all the selected channels. Please only enter the channels numbers for the first wavelength.'};

function vout = nirs_cfg_vout_buildroi1(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});