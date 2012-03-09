function readEprimeOnsets1 = nirs_run_readEprimeOnsets_cfg
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'eprime');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration  Read NIRS onsets from CRIUGM Eprime (Excel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

excelEprime    = cfg_files; %Select raw LOT data files for this subject 
excelEprime.name    = 'Excel from E-Prime'; % The displayed name
excelEprime.tag     = 'excelEprime';       %file names
excelEprime.filter = '.xls';   
excelEprime.num     = [1 Inf];     % Number of inputs required (2D-array with exactly one row and one column)
excelEprime.help    = {'Select excel files imported from E-Prime .edat for the subject.'}; % help text displayed

readEprimeOnsets1      = cfg_exbranch;      
readEprimeOnsets1.name = 'Get E-Prime Onsets';      
readEprimeOnsets1.tag  = 'readEprimeOnsets1'; 
readEprimeOnsets1.val  = {NIRSmat redo1 NIRSmatCopyChoice excelEprime}; 
readEprimeOnsets1.prog = @nirs_run_readEprimeOnsets; 
readEprimeOnsets1.vout = @nirs_cfg_vout_readEprimeOnsets; 
readEprimeOnsets1.help = {'Write over given onset files, permuting onsets as, ',...
    'required so that onsets are in the same order as for the first file. ',...
    'This module should be run by itself, not as part of a larger batch.'};

function vout = nirs_cfg_vout_readEprimeOnsets(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});