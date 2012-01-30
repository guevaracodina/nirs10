function [NIRSmat NIRSmatCopyChoice] = get_common_NIRSmat(NIRSmatReq,dir_name)
NIRSmat         = cfg_files; %Select NIRS.mat for this subject
NIRSmat.name    = 'NIRS.mat'; % The displayed name
NIRSmat.tag     = 'NIRSmat';       %file names
NIRSmat.filter  = 'mat';
NIRSmat.ufilter = '^NIRS.mat$';
NIRSmat.num     = [NIRSmatReq Inf];     % Number of inputs required
NIRSmat.help    = {'Select NIRS.mat for the subject(s).'}; % help text displayed

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