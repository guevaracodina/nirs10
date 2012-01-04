function testreconstruct1 = nirs_run_testreconstruct_cfg


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
% test relevance of the reconstructions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

head_shadow         = cfg_files;
head_shadow.name    = 'Select segmented image';
head_shadow.tag     = 'head_shadow';
head_shadow.ufilter = {'.nii'};
head_shadow.num     = [1 Inf];
head_shadow.help    = {'.'};%The head shadow image has been created by MCsegment and should be located in the directory of the anatomical image

layers_opt        = cfg_menu;
layers_opt.name   = 'Choose the degree of complexity of the phantom';
layers_opt.tag    = 'layers_opt';
layers_opt.labels = {'Homogeneous' '2 layers' '5 layers'};
layers_opt.values = {0,1,2};
layers_opt.val    = {1};
layers_opt.help   = {'-- Homogeneous : phantom homogeneous (like head shadow) with properties of gray matter'
    '-- 2 layers : skin and skull gathered in the first layer then the other layers like grey matter and an inclusion in grey matter'
    '-- 5 layers : and an inclusion in grey matter'}';

inclusion        = cfg_menu;
inclusion.name   = 'Choose the number of inclusions';
inclusion.tag    = 'inclusion';
inclusion.labels = {'0' '1' '2'};
inclusion.values = {0,1,2};
inclusion.val    = {1};
inclusion.help   = {'-- first inclusion in grey matter'
    '-- second inclusion in the skin'}';

% Executable Branch
testreconstruct1      = cfg_exbranch;
testreconstruct1.name = 'Test of reconstructions';
testreconstruct1.tag  = 'testreconstruct1';
testreconstruct1.val  = {NIRSmat DelPreviousData NewDirCopyNIRS head_shadow layers_opt inclusion};
testreconstruct1.prog = @nirs_run_testreconstruct;
testreconstruct1.vout = @nirs_cfg_vout_testreconstruct;
testreconstruct1.help = {'Builds a phantom based on subject 53 (Claudine''s study)'}';

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_testreconstruct(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
