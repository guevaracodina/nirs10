function testreconstruct1 = nirs_run_testreconstruct_cfg
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'testrecon');

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
testreconstruct1.val  = {NIRSmat redo1 NIRSmatCopyChoice head_shadow layers_opt inclusion};
testreconstruct1.prog = @nirs_run_testreconstruct;
testreconstruct1.vout = @nirs_cfg_vout_testreconstruct;
testreconstruct1.help = {'Builds a phantom based on subject 53 (Claudine''s study)'}';

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_testreconstruct(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
