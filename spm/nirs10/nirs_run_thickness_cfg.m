function Uthickness = nirs_run_thickness_cfg

NIRSmat         = cfg_files; %Select NIRS.mat for this subject
NIRSmat.name    = 'NIRS.mat'; % The displayed name
NIRSmat.tag     = 'NIRSmat';       %file names
NIRSmat.filter  = 'mat';
NIRSmat.ufilter = '^NIRS.mat$';
NIRSmat.num     = [1 Inf];     % Number of inputs required
NIRSmat.help    = {'Select NIRS.mat for the subject(s).'}; % help text displayed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Anatomical image utilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SP         = cfg_entry;
SP.name    = 'Point for sampling';
SP.tag     = 'SP';
SP.strtype = 's';
SP.num     = [0 Inf];
SP.val     = {''};
SP.help    = {'Enter the S** or D** or any name of a point registered thanks to Brainsight.'}';

cs_in         = cfg_files;
cs_in.tag     = 'dir_in';
cs_in.name    = 'MonteCarlo output directory';
cs_in.help    = {'Select the MonteCarlo simulation output directory.'};
cs_in.filter = 'dir';
cs_in.val{1} = {''};
cs_in.ufilter = '.*';
cs_in.num     = [0 1];

Uthickness      = cfg_exbranch;
Uthickness.name = 'Thickness below selected points';
Uthickness.tag  = 'Uthickness';
Uthickness.val  = {NIRSmat cs_in SP};
Uthickness.prog = @nirs_run_thickness;
Uthickness.vout = @nirs_cfg_vout_thickness;
Uthickness.help = {'.'};

function vout = nirs_cfg_vout_thickness(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
