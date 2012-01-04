function boxy1 = nirs_run_boxy_cfg

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
%Configuration for BOXY   (boxy1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

input1         = cfg_files; %Select raw BOXY data files for this subject
input1.name    = 'Select BOXY files'; % The displayed name
input1.tag     = 'fnames';       %file names
input1.ufilter = '.0*';    %BOXY files are labeled .001, .002, .003, etc.
%and it is very unlikely that there are more than 99 of them
input1.num     = [1 Inf];     % Number of inputs required
input1.help    = {'Select raw BOXY data files for this subject.'};

age1         = cfg_entry;
age1.name    = 'Subject age';
age1.tag     = 'age1';
age1.strtype = 'r';
age1.num     = [1 1];
age1.def     = @(val)nirs_get_defaults('readNIRS.boxy1.generic1.subj.age1', val{:});
age1.help    = {'Age of the subject. Used later for OD to HbO/HbR conversion.'};

raw_onset_files        = cfg_files;
raw_onset_files.name    = 'Select onset files'; % The displayed name
raw_onset_files.tag     = 'raw_onset_files';
raw_onset_files.num     = [0 Inf];     % Number of inputs required
raw_onset_files.val{1}  = {''};
raw_onset_files.help    = {'Optional: Select raw onset files. '
    'Can be added at a later stage.'
    'Must specify one file for each data file, in same order.'}'; % help text displayed

anatT1         = cfg_files; %Select T1 for this subject
anatT1.name    = 'Raw anatomical image (optional)'; % The displayed name
anatT1.tag     = 'anatT1';       %file names
anatT1.filter  = 'image';
anatT1.ufilter = '.*';
anatT1.val{1}  = {''};
anatT1.num     = [0 Inf];     % Number of inputs required
anatT1.help    = {'CORRECTION: THIS IS NO LONGER OPTIONAL.'
    'It is required for the coregistration module and the contrast module'
    'Optional, can be specified in MC Segment, or earlier '
    'and be available in NIRS.mat structure.'
    'Select raw anatomical image(s) for the subject(s). '
    'If several subjects, the images '
    'must be in the same order as the NIRS.mat structures.'}';

subj         = cfg_branch;
subj.tag     = 'subj';
subj.name    = 'Subject';
subj.val     = {input1 age1 raw_onset_files anatT1};
subj.help    = {'This module allows multi-subject processing, '
    'generating a NIRS.mat file for each subject. '
    'Note that a list of links to the NIRS.mat structures will be '
    'available as a virtual output for further processing'}';

generic1         = cfg_repeat;
generic1.tag     = 'generic1';
generic1.name    = 'Subject';
generic1.help    = {'Note simple data organization that is recommended: '
    '1- directories with files for each subject must be in the same root directory. '
    '2- .prj file must have same name as expname and in folder \mtg\ '
    '3- anatomic files must be in folder \T1\'}';
generic1.values  = {subj};
generic1.num     = [1 Inf];

prj_path         = cfg_entry; %path
prj_path.def     = @(val)nirs_get_defaults('readNIRS.boxy1.config_path.prj_path', val{:});
prj_path.name    = 'path for .prj file';
prj_path.tag     = 'prj_path';
prj_path.strtype = 's';
prj_path.num     = [1 Inf];
%prj_path.val     = {'mtg'}; %Not using .def due to a bug, which is not
%understood: the default value is not found
prj_path.help    = {'Path for .prj file: should be something like ..\mtg\ (omit backslashes). '
    'Hint: if default entry (mtg) is missing, close spm and restart (without '
    'closing Matlab or clearing variables in the workspace.'}';

%path structure
T1_path         = cfg_entry; %path
T1_path.name    = 'path for anatomical files';
T1_path.tag     = 'T1_path';
T1_path.strtype = 's';
T1_path.num     = [1 Inf];
T1_path.def     = @(val)nirs_get_defaults('readNIRS.boxy1.config_path.T1_path', val{:});
T1_path.help    = {'Path for T1 file: should be something like ..\T1\ (omit backslashes)'};

%path structure
output_path         = cfg_entry; %path
output_path.name    = 'path for .nir output files';
output_path.tag     = 'output_path';
output_path.strtype = 's';
output_path.num     = [1 Inf];
output_path.def     = @(val)nirs_get_defaults('readNIRS.boxy1.config_path.output_path', val{:});
output_path.help    = {'Path for .nir output files: should be something like ..\dataSPM\ (omit backslashes)'};

config_path         = cfg_branch;
config_path.tag     = 'config_path';
config_path.name    = 'Path Configuration options';
config_path.val     = {prj_path T1_path output_path};
config_path.help    = {''};

%Light wavelengths
Lambda         = cfg_entry; %Lambda
Lambda.name    = 'Laser wavelengths';
Lambda.tag     = 'Lambda';
Lambda.strtype = 'r';
Lambda.num     = [1 Inf];
Lambda.def     = @(val)nirs_get_defaults('readNIRS.boxy1.cf1.Lambda', val{:});
Lambda.help    = {'Near Infrared laser wavelengths. Note order is critical and'
    'must correspond to recording order in raw data files.'}';

%wavelengths sensitive to HbO - NOT USED - calculated in the code based on
%info stored in NIRS.Cf.dev.wl
% LambdaHbO         = cfg_entry; %Lambda
% LambdaHbO.name    = 'Wavelengths most sensitive to HbO';
% LambdaHbO.tag     = 'LambdaHbO';
% LambdaHbO.strtype = 'r';
% LambdaHbO.num     = [1 Inf];
% LambdaHbO.def     = @(val)nirs_get_defaults('readNIRS.boxy1.cf1.LambdaHbO', val{:});
% LambdaHbO.help    = {'Enter a vector, with one entry for each wavelength,'
%             'Indicating with a 1 if the wavelength is most sensitive to HbO'
%             'and with a 0 if it is less sensitive. With more than 2 wavelengths,'
%             'several wavelengths could be sensitive to HbO. This vector of '
%             'Booleans will be used for detection of heart rate and Mayer waves,'
%             'and for the configuration of Monte Carlo files.'}';

%Input frequency
input2 = cfg_entry;
input2.name    = 'Input frequency';
input2.tag     = 'freq';
input2.strtype = 'r';
input2.num     = [1 1];     % Number of inputs required (2D-array with exactly one row and one column)
input2.def     = @(val)nirs_get_defaults('readNIRS.boxy1.cf1.freq', val{:});
input2.help    = {'Input data frequency in Hertz.'};

%Minimum distance
distmin         = cfg_entry;
distmin.name    = 'Minimum distance';
distmin.tag     = 'distmin';
distmin.strtype = 'r';
distmin.num     = [1 1];
distmin.def     = @(val)nirs_get_defaults('readNIRS.boxy1.cf1.distmin', val{:});
distmin.help    = {'Cutoff: Minimum Cartesian channel distance in centimeters.'};

%Maximum distance
distmax         = cfg_entry;
distmax.name    = 'Maximum distance';
distmax.tag     = 'distmax';
distmax.strtype = 'r';
distmax.num     = [1 1];
distmax.def     = @(val)nirs_get_defaults('readNIRS.boxy1.cf1.distmax', val{:});
distmax.help    = {'Cutoff: Maximum Cartesian channel distance in centimeters.'};

%Size of blocks used for processing
sizebloc         = cfg_entry;
sizebloc.name    = 'Block size';
sizebloc.tag     = 'sizebloc';
sizebloc.strtype = 'r';
sizebloc.num     = [1 1];
sizebloc.def     = @(val)nirs_get_defaults('readNIRS.boxy1.cf1.sizebloc', val{:});
sizebloc.help    = {'Size of blocks used for processing whole file in chunks.'};

%Number of MUX
nb_Mux         = cfg_entry; %nb_Mux
nb_Mux.name    = 'MUX number';
nb_Mux.tag     = 'nb_Mux';
nb_Mux.strtype = 'r';
nb_Mux.num     = [1 1];
nb_Mux.def     = @(val)nirs_get_defaults('readNIRS.boxy1.cf1.nb_Mux', val{:});
nb_Mux.help    = {'Number of MUX (Note: NOT the same as the number of sources!)'};

%MaxSources
MaxSources         = cfg_entry; %MaxSources
MaxSources.name    = 'Maximum number of sources';
MaxSources.tag     = 'MaxSources';
MaxSources.strtype = 'r';
MaxSources.num     = [1 1];
MaxSources.def     = @(val)nirs_get_defaults('readNIRS.boxy1.cf1.MaxSources', val{:});
MaxSources.help    = {'Maximum number of sources'};

%Maximum number of detectors
nb_Det         = cfg_entry; %nb_Det
nb_Det.name    = 'Maximum number of detectors';
nb_Det.tag     = 'nb_Det';
nb_Det.strtype = 'r';
nb_Det.num     = [1 1];
nb_Det.def     = @(val)nirs_get_defaults('readNIRS.boxy1.cf1.nb_Det', val{:});
nb_Det.help    = {'Maximum number of detectors'};

%Maximum number of detectors
MaxElectrodes         = cfg_entry; %MaxElectrodes
MaxElectrodes.name    = 'Maximum number of electrodes';
MaxElectrodes.tag     = 'MaxElectrodes';
MaxElectrodes.strtype = 'r';
MaxElectrodes.num     = [1 1];
MaxElectrodes.def     = @(val)nirs_get_defaults('readNIRS.boxy1.cf1.MaxElectrodes', val{:});
MaxElectrodes.help    = {'Maximum and habitual number of electrodes. This is only used'
    'in the interpolation, when one needs to fill in for missing electrode positions'}';

%use10_10system
use10_10system      = cfg_menu;
use10_10system.tag  = 'use10_10system';
use10_10system.name = 'Use 10 10 system';
use10_10system.labels = {'True','False: use 10 20 system'};
use10_10system.values = {1,0};
use10_10system.def  = @(val)nirs_get_defaults('readNIRS.boxy1.cf1.use10_10system', val{:});
use10_10system.help = {'10 10 system allows more precise location of optodes'
    'for viewing in .nir files.' }';

%Downsampling factor
resample         = cfg_entry; %resample
resample.name    = 'Downsampling factor';
resample.tag     = 'resample';
resample.strtype = 'r';
resample.num     = [1 1];
resample.def     = @(val)nirs_get_defaults('readNIRS.boxy1.cf1.resample', val{:});
resample.help    = {'Downsampling factor: use 1 to keep all data.'
    'Note that no filtering will be done at this stage;'
    'therefore beware of aliasing artefacts.'}';

save_bin1      = cfg_menu;
save_bin1.tag  = 'save_bin1';
save_bin1.name = 'Save binary files';
save_bin1.labels = {'True','False'};
save_bin1.values = {1,0};
save_bin1.def  = @(val)nirs_get_defaults('readNIRS.boxy1.cf1.save_bin1', val{:});
save_bin1.help = {'Write out binary files for NIRS data.'
    'Note that binary files are much faster to access than ASCII or .mat files'}';

cf1         = cfg_branch;
cf1.tag     = 'cf1';
cf1.name    = 'Configuration options';
cf1.val     = {Lambda input2 distmin distmax save_bin1 sizebloc ...
    nb_Mux MaxSources nb_Det MaxElectrodes use10_10system resample};
cf1.help    = {'Configuration used for all subjects to be preprocessed.'};

% Executable Branch
boxy1      = cfg_exbranch;
boxy1.name = 'ReadBoxy';
boxy1.tag  = 'boxy1';
boxy1.val  = {generic1 config_path cf1};
boxy1.prog = @nirs_run_boxy;
boxy1.vout = @nirs_cfg_vout_boxy;
boxy1.help = {'Select raw BOXY data files for this subject.'};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_boxy(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
