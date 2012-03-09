function [NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(NIRSmatReq,dir_name)
redo1      = cfg_menu;
redo1.tag  = 'force_redo';
redo1.name = 'Force processing';
redo1.labels = {'False','True'};
redo1.values = {0,1};
redo1.val  = {0};
redo1.help = {'Force redoing this processing even when it has been done already'};

NIRSmat         = cfg_files; %Select NIRS.mat for this subject
if NIRSmatReq
    NIRSmat.name    = 'NIRS.mat'; % The displayed name
else
    NIRSmat.name    = 'NIRS.mat (optional)'; % The displayed name    
end
NIRSmat.help    = {'Select NIRS.mat for the subject(s).'}; % help text displayed
NIRSmat.tag     = 'NIRSmat';       %file names
NIRSmat.filter  = 'mat';
NIRSmat.ufilter = '^NIRS.mat$';
NIRSmat.num     = [NIRSmatReq Inf];     % Number of inputs required

NIRSmatOverwrite         = cfg_branch;
NIRSmatOverwrite.tag     = 'NIRSmatOverwrite';
NIRSmatOverwrite.name    = 'Overwrite NIRS.mat structure'; 
NIRSmatOverwrite.help    = {'Will not copy NIRS structure.'
            'This will write over the previous NIRS.mat'}';

NewNIRSdir         = cfg_entry;
NewNIRSdir.name    = 'New directory for NIRS.mat';
NewNIRSdir.tag     = 'NewNIRSdir';       
NewNIRSdir.strtype = 's';
NewNIRSdir.val{1}    = dir_name;
NewNIRSdir.num     = [1 Inf];     
NewNIRSdir.help    = {'Directory for NIRS.mat.'}'; 

NIRSmatCopy         = cfg_branch;
NIRSmatCopy.tag     = 'NIRSmatCopy';
NIRSmatCopy.name    = 'Create new directory and copy NIRS structure there'; 
NIRSmatCopy.val     = {NewNIRSdir};
NIRSmatCopy.help    = {'Create new directory and copy NIRS structure there.'}';
        
%Common to most modules: for creating a new directory and copying NIRS.mat
NIRSmatCopyChoice           = cfg_choice;
NIRSmatCopyChoice.name      = 'Choose NIRS copy method';
NIRSmatCopyChoice.tag       = 'NIRSmatCopyChoice';
NIRSmatCopyChoice.values    = {NIRSmatOverwrite NIRSmatCopy}; 
NIRSmatCopyChoice.val       = {NIRSmatOverwrite}; 
NIRSmatCopyChoice.help      = {'Choose whether to overwrite the NIRS.mat structure'
            'or to create a new directory'
            'and copy the NIRS.mat structure there'}'; 