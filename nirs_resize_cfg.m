function resize1 = nirs_resize_cfg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

image_in         = cfg_files;
image_in.name    = 'Anatomical segmented image';
image_in.tag     = 'image_in';
image_in.filter  = 'image';
image_in.num     = [1 1];
image_in.help    = {'Select the SAME image as the one used to run the MC simulations'};

out_dir         = cfg_files;
out_dir.tag     = 'out_dir';
out_dir.name    = 'Output Directory';
out_dir.help    = {'Select a directory where the output image will be written.'};
out_dir.filter = 'dir';
out_dir.ufilter = '.*';
out_dir.num     = [1 1];

out_dim      = cfg_entry;
out_dim.tag  = 'out_dim';
out_dim.name = 'Output dimension';
out_dim.val = {1};
out_dim.strtype = 'r';
out_dim.num     = [1 3];
out_dim.def  = @(val)nirs_get_defaults('coregNIRS.resize1.out_dim', val{:});
out_dim.help = {['Enter output image size or let [1 1 1] if you just ',...
    'want to get isotropic voxels image.']};

out_dt      = cfg_entry;
out_dt.tag  = 'out_dt';
out_dt.name = 'Data type';
out_dt.val = {1};
out_dt.strtype = 's';
out_dt.num     = [1 Inf];
out_dt.def  = @(val)nirs_get_defaults('coregNIRS.resize1.out_dt', val{:});
out_dt.help = {['Enter output image size or let ''same'' if you just ',...
    'want to get isotropic voxels image.']};

out_autonaming      = cfg_menu;
out_autonaming.tag  = 'out_autonaming';
out_autonaming.name = 'Automatic output naming';
out_autonaming.labels = {'Yes','No'};
out_autonaming.values = {0,1};
out_autonaming.def  = @(val)nirs_get_defaults('coregNIRS.resize1.out_autonaming', val{:});
out_autonaming.help = {['Choose wheather you want to choose the name ',...
    'of the output or not. If answer is ''Yes'', please change enter name.']};

out_prefix      = cfg_entry;
out_prefix.tag  = 'out_prefix';
out_prefix.name = 'Prefix of the output image (resize1)';
out_prefix.strtype = 's';
out_prefix.num     = [1 Inf];
out_prefix.help = {['You can choose to give a particular prefix to the ',...
    'output image. This prefix will be added at the left of the name ',...
    'of the image. A default name will be given otherwise.']};

% Executable Branch
resize1      = cfg_exbranch;
resize1.name = 'Resize image';
resize1.tag  = 'resize1';
resize1.val  = {image_in out_dir out_dim out_dt out_autonaming out_prefix};
resize1.prog = @nirs_resize;
resize1.vout = @nirs_cfg_vout_resize;
resize1.help = {'Resize the input image with respect to output size.'};


function vout = nirs_cfg_vout_resize(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});