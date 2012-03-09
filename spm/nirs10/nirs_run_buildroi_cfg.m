function buildroi1 = nirs_run_buildroi_cfg
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'ROI');

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
buildroi1.val  = {NIRSmat redo1 NIRSmatCopyChoice keepAllChannels image_in output_prefix};
buildroi1.prog = @nirs_run_buildroi2;
buildroi1.vout = @nirs_cfg_vout_buildroi1;
buildroi1.help = {'Define region of interest containing all the selected channels. Please only enter the channels numbers for the first wavelength.'};

function vout = nirs_cfg_vout_buildroi1(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});