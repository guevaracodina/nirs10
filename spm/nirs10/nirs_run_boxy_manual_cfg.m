function boxy_manual1 = nirs_run_boxy_manual_cfg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration for BOXY MANUAL  (boxy1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Apath        = cfg_files;
Apath.name    = 'Analysis path';
Apath.tag     = 'Apath';       %file names
Apath.filter  = 'dir';
Apath.ufilter = '.*';
Apath.num     = [1 1];
Apath.help    = {'Choose folder where you want to put the analysis.'
    'If the folder does not exist yet, you must first create it '
    'before you can select it.'}';

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

raw_onset_files2        = cfg_files;
raw_onset_files2.name    = 'Select onset files'; % The displayed name
raw_onset_files2.tag     = 'raw_onset_files2';
raw_onset_files2.num     = [1 Inf];     % Number of inputs required
raw_onset_files2.val{1}  = {''};
raw_onset_files2.help    = {'Select raw onset files. '
    'Must specify one file for each data file, in same order.'}'; % help text displayed

anatT1         = cfg_files; %Select T1 for this subject
anatT1.name    = 'Raw anatomical image (optional)'; % The displayed name
anatT1.tag     = 'anatT1';       %file names
anatT1.filter  = 'image';
anatT1.ufilter = '.*';
anatT1.val{1}  = {''};
anatT1.num     = [0 Inf];     % Number of inputs required
anatT1.help    = {'Optional, can be specified in MC Segment, or earlier '
    'and be available in NIRS.mat structure.'
    'Select raw anatomical image(s) for the subject(s). '
    'If several subjects, the images '
    'must be in the same order as the NIRS.mat structures.'}';

prjFile        = cfg_files;
prjFile.name    = 'Montage file'; % The displayed name
prjFile.tag     = 'prjFile';       %file names
prjFile.ufilter = '.prj';
prjFile.num     = [1 1];     % Number of inputs required
prjFile.help    = {'Select montage file (.prj), from MTG folder.'}';

subj2         = cfg_branch;
subj2.tag     = 'subj2';
subj2.name    = 'Subject';
subj2.val     = {Apath input1 prjFile age1 raw_onset_files2 anatT1};
subj2.help    = {'This module allows multi-subject processing, '
    'generating a NIRS.mat file for each subject. '
    'Note that a list of links to the NIRS.mat structures will be '
    'available as a virtual output for further processing'}';

generic2         = cfg_repeat;
generic2.tag     = 'generic2';
generic2.name    = 'Subject';
generic2.help    = {'For this module, more flexibility in location of files'
    'is allowed, compared to the previous Read BOXY module'}';
generic2.values  = {subj2};
generic2.num     = [1 Inf];

% %path structure
% T1_path         = cfg_entry; %path
% T1_path.name    = 'path for anatomical files';
% T1_path.tag     = 'T1_path';
% T1_path.strtype = 's';
% T1_path.num     = [1 Inf];
% T1_path.def     = @(val)nirs_get_defaults('readNIRS.boxy1.config_path.T1_path', val{:});
% T1_path.help    = {'Path for T1 file: should be something like ..\T1\ (omit backslashes)'};
%
% %path structure
% output_path         = cfg_entry; %path
% output_path.name    = 'path for .nir output files';
% output_path.tag     = 'output_path';
% output_path.strtype = 's';
% output_path.num     = [1 Inf];
% output_path.def     = @(val)nirs_get_defaults('readNIRS.boxy1.config_path.output_path', val{:});
% output_path.help    = {'Path for .nir output files: should be something like ..\dataSPM\ (omit backslashes)'};

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

config_path2         = cfg_branch;
config_path2.tag     = 'config_path2';
config_path2.name    = 'Path Configuration options';
config_path2.val     = {T1_path output_path};
config_path2.help    = {''};

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
boxy_manual1      = cfg_exbranch;
boxy_manual1.name = 'ReadBoxyManual';
boxy_manual1.tag  = 'boxy_manual1';
boxy_manual1.val  = {generic2 config_path2 cf1};
boxy_manual1.prog = @nirs_run_boxy_manual;
boxy_manual1.vout = @nirs_cfg_vout_boxy_manual;
boxy_manual1.help = {'Select raw BOXY data files for this subject.'};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_boxy_manual(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});


