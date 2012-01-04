function permuteOnsets1 = nirs_run_permuteOnsets_cfg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration  Permute Onsets if required
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Utility to permute onsets so that
onset_files        = cfg_files;
onset_files.name    = 'Select onset files to permute';
onset_files.tag     = 'onset_files';
onset_files.filter  = 'mat';
onset_files.num     = [1 Inf];
onset_files.help    = {'Select onset files to be permuted to match order of onsets of first file.'}; % help text displayed

% Executable Branch
permuteOnsets1      = cfg_exbranch;
permuteOnsets1.name = 'Permute Onsets';
permuteOnsets1.tag  = 'permuteOnsets1';
permuteOnsets1.val  = {onset_files};
permuteOnsets1.prog = @nirs_run_permuteOnsets;
permuteOnsets1.vout = @nirs_cfg_vout_permuteOnsets;
permuteOnsets1.help = {'Write over given onset files, permuting onsets as, ',...
    'required so that onsets are in the same order as for the first file. ',...
    'This module should be run by itself, not as part of a larger batch.'};

function vout = nirs_cfg_vout_permuteOnsets(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
