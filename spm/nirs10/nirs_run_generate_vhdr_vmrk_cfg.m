function generate_vhdr_vmrk = nirs_run_generate_vhdr_vmrk_cfg

NIRSmat         = cfg_files; %Select NIRS.mat for this subject
NIRSmat.name    = 'NIRS.mat'; % The displayed name
NIRSmat.tag     = 'NIRSmat';       %file names
NIRSmat.filter  = 'mat';
NIRSmat.ufilter = '^NIRS.mat$';
NIRSmat.num     = [1 Inf];     % Number of inputs required
NIRSmat.help    = {'Select NIRS.mat for the subject(s).'}; % help text displayed


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4.7 Generate header and marker files for Analyzer based on NIRS.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Executable Branch
generate_vhdr_vmrk      = cfg_exbranch;
generate_vhdr_vmrk.name = 'Generate header and marker files';
generate_vhdr_vmrk.tag  = 'generate_vhdr_vmrk';
generate_vhdr_vmrk.val  = {NIRSmat};
generate_vhdr_vmrk.prog = @nirs_run_generate_vhdr_vmrk;
generate_vhdr_vmrk.vout = @nirs_cfg_vout_generate_vhdr_vmrk;
generate_vhdr_vmrk.help = {'Generate Brain Vision Analyzer ',...
    'header (.vhrd) and marker (.vmrk) files based on NIRS structure.'};

function vout = nirs_cfg_vout_generate_vhdr_vmrk(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});