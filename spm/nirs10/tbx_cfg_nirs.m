function nirs10 = tbx_cfg_nirs
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________

addpath(fileparts(which(mfilename)));

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
age1.def     = @(val)nirs_get_defaults('nirs10.readNIRS.boxy1.generic.subj.age1', val{:}); 
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
anatT1.help    = {'Optional, can be specified in MC Segment, or earlier '
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

generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Subject';
generic.help    = {'Note simple data organization that is recommended: '
    '1- directories with files for each subject must be in the same root directory. '
    '2- .prj file must have same name as expname and in folder \mtg\ '
    '3- anatomic files must be in folder \T1\'}';
generic.values  = {subj};
generic.num     = [1 Inf];

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
input2         = cfg_entry; 
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
boxy1.val  = {generic config_path cf1};   
boxy1.prog = @nirs_run_boxy;  
boxy1.vout = @nirs_cfg_vout_boxy;
boxy1.help = {'Select raw BOXY data files for this subject.'};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_boxy(job)
vout = cfg_dep;                    
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration for IUGM (Techen CW5 [UNF] or CW6 [LESCA])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TopoData        = cfg_files;
TopoData.tag    = 'TopoData';
TopoData.name   = 'TopoData';
TopoData.filter = '.mat';
TopoData.num    = [0 1];
TopoData.val{1} = {''};
TopoData.help   = {'Select TopoData if you want to run GLM. You can choose either the template or the one corresponding to your subject if it has already been calculated'};

age1         = cfg_entry;
age1.name    = 'Subject age';
age1.tag     = 'age1';
age1.strtype = 'r';
age1.num     = [1 1];
age1.def     = @(val)nirs_get_defaults('nirs10.readNIRS.criugm1.generic.subj.age1', val{:});
age1.help    = {'Age of the subject. Used later for OD to HbO/HbR conversion.'};

text_brainsight         = cfg_files;
text_brainsight.tag     = 'text_brainsight';
text_brainsight.name    = 'Text file from Brainsight';
text_brainsight.filter  = '.txt';
text_brainsight.ufilter = '.*';
text_brainsight.num     = [1 1];
text_brainsight.help    = {'Select the text file from Brainsight.'};
 
T1_vitamins           = cfg_branch; 
T1_vitamins.name      = 'Vitamins markers on T1';
T1_vitamins.tag       = 'T1_vitamins';
T1_vitamins.help      = {'The helmet will be read in future module from T1 image and positions of vitamins in the image.'};

no_helmet = cfg_branch; 
no_helmet.name      = 'No helmet information';
no_helmet.tag       = 'no_helmet';
no_helmet.help      = {'No helmet information available. Some modules won''t work properly.'};

helmet         = cfg_choice;
helmet.tag     = 'helmet';
helmet.name    = 'Helmet';
helmet.values = {text_brainsight T1_vitamins no_helmet};
helmet.val     = {text_brainsight};
helmet.help    = {'If you choose a Brainsight text file, it will be used to determine all you need about sources, detectors, channels,... Otherwise, information will be extrqcted from the first file.'};

nirs_files         = cfg_files;
nirs_files.name    = '''^.nirs'' files'; % The displayed name
nirs_files.tag     = 'nirs_files';       %file names
nirs_files.filter  = 'nirs';
nirs_files.num     = [0 Inf];     % Number of inputs required
nirs_files.help    = {'Select all the sessions sharing the same device and helmet.'}; % help text displayed

protocol        = cfg_files;
protocol.tag    = 'protocol';
protocol.name   = 'Experimental conditions';
protocol.filter = '.mat';
protocol.num    = [0 Inf];
protocol.val{1} = {''};
protocol.help   = {['Select the "SPM conditions" .mat files. The order must correspond to that of the .nirs files.' ...
    ' If some (N) .nirs file have no corresponding protocol, they must be selected as the N last ones.' ...
    ' The format of the "SPM conditions" files is as follows: for n conditions, 3 (1xn) cell arrays contain' ...
    ' the names (string), onsets (vector) and durations (1 scalar or vector), ...' ...
    ' (where length(onsets{i}) = number of events for condition i.']};

CWsystem      = cfg_menu;
CWsystem.tag  = 'CWsystem';
CWsystem.name = 'CW system used';
CWsystem.labels = {'CW5','CW6'};
CWsystem.values = {5,6};
CWsystem.def  = @(val)nirs_get_defaults('readNIRS.criugm1.CWsystem', val{:});
CWsystem.help = {'Help'};

% dataset      = cfg_branch;
% dataset.tag  = 'dataset';
% dataset.name = 'Data set';
% dataset.val  = {nirs_file helmet CWsystem};
% dataset.help = {'A data set is a set whose ''.nirs'' files (NIRS raw data files) shares the same Helmet and CW device.'};

% NIRSsessions         = cfg_repeat;
% NIRSsessions.tag     = 'NIRSsessions';
% NIRSsessions.name    = 'NIRS sessions';
% NIRSsessions.help    = {'Help'};
% NIRSsessions.values  = {nirs_files};
% NIRSsessions.num     = [1 Inf];
% NIRSsessions.help = {'Enter all the data you want to process for the selected subject. These data must be sorted by helmet and device. Create as many sets as the number of different configurations.'};

baseline_method = cfg_menu;
baseline_method.tag  = 'baseline_method';
baseline_method.name = 'Method to calculate Baseline';
baseline_method.labels = {'Median','Initial Value','Mean'};
baseline_method.values = {0,1,2};
baseline_method.def  = @(val)nirs_get_defaults('readNIRS.criugm1.baseline_method', val{:});
baseline_method.help = {'Help'};

subj_path         = cfg_files;
subj_path.tag     = 'subj_path';
subj_path.name    = 'Subject Directory';
subj_path.help    = {'Select a directory where the NIRS matrix will be written.'};
subj_path.filter = 'dir';
subj_path.ufilter = '.*';
subj_path.num     = [1 1];

subj         = cfg_branch;
subj.tag     = 'subj';
subj.name    = 'Subject';
subj.val     = {age1 subj_path anatT1 helmet TopoData CWsystem protocol nirs_files baseline_method};
subj.help    = {'Subject'};

generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Subjects';
generic.help    = {'Help'};
generic.values  = {subj};
generic.num     = [1 Inf];

% The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
% Executable Branch
criugm1      = cfg_exbranch;
criugm1.name = 'Read and format CRIUGM data';
criugm1.tag  = 'criugm1';
criugm1.val  = {generic};
criugm1.prog = @nirs_run_criugm;
criugm1.vout = @nirs_cfg_vout_criugm;
criugm1.help = {'Help'};

%make NIRS.mat available as a dependency
    function vout = nirs_cfg_vout_criugm(job)
        vout = cfg_dep;                     % The dependency object
        vout.sname      = 'NIRS.mat';
        vout.src_output = substruct('.','NIRSmat');
        vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration for LOT   (lot1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

input1         = cfg_files; %Select raw LOT data files for this subject 
input1.name    = 'Select LOT files'; % The displayed name
input1.tag     = 'fnames';       %file names
%input1.ufilter = '.0*';   
input1.num     = [1 Inf];     % Number of inputs required (2D-array with exactly one row and one column)
input1.help    = {'Select raw LOT data files for this subject.'}; % help text displayed

% Executable Branch
lot1      = cfg_exbranch;       
lot1.name = 'ReadLot';             
lot1.tag  = 'lot1'; 
lot1.val  = {input1};   
lot1.prog = @nirs_run_lot;  
lot1.vout = @nirs_cfg_vout_lot; 
lot1.help = {'Select raw LOT data files for this subject.'};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_lot(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration  Read NIRS onsets to generate input to General Linear Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration  Read NIRS onsets for epilepsy from Analyzer 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq_NIRS1      = cfg_entry;
freq_NIRS1.tag  = 'freq_NIRS1';
freq_NIRS1.name = 'Frequency of NIRS data for GLM';
freq_NIRS1.val{1} = []; %{1.9531}; %{19.5312};
freq_NIRS1.strtype = 'r';  
freq_NIRS1.num     = [0 Inf]; 
freq_NIRS1.help    = {'Specify frequency of NIRS data (optional).'}; 

dp_NIRS1      = cfg_entry;
dp_NIRS1.tag  = 'dp_NIRS1';
dp_NIRS1.name = 'Number of data points in NIRS data for GLM';
dp_NIRS1.val{1} = []; %{1758}; 
dp_NIRS1.strtype = 'r';  
dp_NIRS1.num     = [0 Inf]; 
dp_NIRS1.help    = {'Specify number of data time points in NIRS data for the GLM (optional).'}; 

NIRSmat_optional         = cfg_files; %Select NIRS.mat for this subject 
NIRSmat_optional.name    = 'NIRS.mat'; % The displayed name
NIRSmat_optional.tag     = 'NIRSmat_optional';       %file names
NIRSmat_optional.filter  = 'mat';
NIRSmat_optional.ufilter = '^NIRS.mat$';  
NIRSmat_optional.val{1}  = {''};
NIRSmat_optional.num     = [0 Inf];     % Number of inputs required 
NIRSmat_optional.help    = {'Select NIRS.mat for the subject(s).'
    'If selecting more than one NIRS.mat, the onsets must have been already '
    'specified at an earlier stage. Otherwise, only one subject can be processed.'}'; % help text displayed


% Executable Branch
AnalyzerOnsets      = cfg_exbranch;       
AnalyzerOnsets.name = 'Read NIRS onsets';            
AnalyzerOnsets.tag  = 'AnalyzerOnsets'; 
AnalyzerOnsets.val  = {NIRSmat_optional raw_onset_files freq_NIRS1 dp_NIRS1}; 
AnalyzerOnsets.prog = @nirs_run_AnalyzerOnsets;  
AnalyzerOnsets.vout = @nirs_cfg_vout_AnalyzerOnsets; 
AnalyzerOnsets.help = {'Select NIRS structures (optional) and/or '
    'Analyzer 2 export files of onsets (also optional), '
    'to generate files of onsets and of confound regressors (pulse). '
    'This module can now be run by itself or as part of a larger batch.'
    'If specifying onset files, only one subject should be run.'}';

function vout = nirs_cfg_vout_AnalyzerOnsets(job)
vout = cfg_dep;                    
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration  Read NIRS onsets from CRIUGM Eprime (Excel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ep_output         = cfg_files;
ep_output.name    = 'Eprime output (Excel file)';
ep_output.tag     = 'ep_output';
ep_output.num     = [1 Inf];
ep_output.help    = {'Select Excel file output from Eprime.'}';

columns         = cfg_entry;
columns.tag     = 'columns';
columns.name    = 'Columns to read';
columns.strtype = 's';
columns.num     = [1 Inf];
columns.def     = @(val)nirs_get_defaults('readOnsets.readEprimeOnsets.columns', val{:});
columns.help = {''};

% Executable Branch
readEprimeOnsets      = cfg_exbranch;       
readEprimeOnsets.name = 'Read Eprime onsets';            
readEprimeOnsets.tag  = 'readEprimeOnsets'; 
readEprimeOnsets.val  = {ep_output columns}; 
readEprimeOnsets.prog = @nirs_run_readEprimeOnsets;  
readEprimeOnsets.vout = @nirs_cfg_vout_readEprimeOnsets; 
readEprimeOnsets.help = {'Read Eprime output data.'}';

function vout = nirs_cfg_vout_readEprimeOnsets(job)
vout = cfg_dep;                    
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

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
permuteOnsets      = cfg_exbranch;      
permuteOnsets.name = 'Permute Onsets';      
permuteOnsets.tag  = 'permuteOnsets'; 
permuteOnsets.val  = {onset_files}; 
permuteOnsets.prog = @nirs_run_permuteOnsets; 
permuteOnsets.vout = @nirs_cfg_vout_permuteOnsets; 
permuteOnsets.help = {'Write over given onset files, permuting onsets as, ',...
    'required so that onsets are in the same order as for the first file. ',...
    'This module should be run by itself, not as part of a larger batch.'};

function vout = nirs_cfg_vout_permuteOnsets(job)
vout = cfg_dep;                     
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tests: add stimuli and corresponding HRFs to raw NIRS signals   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

testDupChannels         = cfg_menu;
testDupChannels.tag     = 'testDupChannels';
testDupChannels.name    = 'Duplicate number of channels';
testDupChannels.help    = {'For each specified channel, add stimuli.'
        'and also add a copy of that channel, without stimuli.'}';
testDupChannels.labels = {
               'True'
               'False'
                }';
testDupChannels.values = {1, 0};
testDupChannels.def     = @(val)nirs_get_defaults(...
    'readOnsets.addTestStimuli.testDupChannels', val{:}); 

testStimulusName         = cfg_entry;
testStimulusName.name    = 'Name of the stimulus';
testStimulusName.tag     = 'testStimulusName';
testStimulusName.strtype = 's';
testStimulusName.num     = [1 Inf];
testStimulusName.def  = @(val)nirs_get_defaults(...
    'readOnsets.addTestStimuli.testStimulusName', val{:});
testStimulusName.help    = {'Name of the stimulus.'};

keepChannels         = cfg_entry; 
keepChannels.name    = 'List of channels to keep';
keepChannels.tag     = 'keepChannels';       
keepChannels.strtype = 'r';
keepChannels.num     = [1 Inf];     
keepChannels.def     = @(val)nirs_get_defaults('readOnsets.addTestStimuli.keepChannels', val{:}); 
keepChannels.help    = {'Enter channel numbers to keep.'}; 

AllChannels           = cfg_branch; 
AllChannels.name      = 'Keep ALL channels';
AllChannels.tag       = 'AllChannels';
AllChannels.help      = {'Keep all channels'}; 
   
keepAllChannels        = cfg_choice;
keepAllChannels.name   = 'Specify channels to keep';
keepAllChannels.tag    = 'keepAllChannels';
keepAllChannels.values = {AllChannels keepChannels};
%Do not know how to specify a default value by a call using .def for a 
%cfg_choice object
keepAllChannels.val    = {keepChannels};
keepAllChannels.help   = {'Choose whether to keep all channels or select a subset.'};

%Which channels to test
testChannels         = cfg_entry; 
testChannels.name    = 'Test channel numbers';
testChannels.tag     = 'testChannels';       
testChannels.strtype = 'r';
testChannels.num     = [1 Inf];     
testChannels.def     = @(val)nirs_get_defaults('readOnsets.addTestStimuli.testChannels', val{:}); 
testChannels.help    = {'Enter channel numbers where the stimuli will be added.'
    'It is sufficient to give only the channel numbers for the first '
    'wavelength, and the code will add the mirror channels for the other wavelengths'}'; 

%Test stimuli number
testStimuliNumber         = cfg_entry; 
testStimuliNumber.name    = 'Test stimuli number';
testStimuliNumber.tag     = 'testStimuliNumber';       
testStimuliNumber.strtype = 'r';
testStimuliNumber.num     = [1 1];     
testStimuliNumber.def     = @(val)nirs_get_defaults(...
    'readOnsets.addTestStimuli.testStimuliNumber', val{:}); 
testStimuliNumber.help    = {'Enter number of stimuli to be added.'}; 

testBP           = cfg_branch; %empty branch
testBP.name      = 'Block Paradigm';
testBP.tag       = 'testBP';
testBP.val       = {};
testBP.help      = {'Block Paradigm'}; 

testSeed1         = cfg_entry; 
testSeed1.name    = 'Random seed';
testSeed1.tag     = 'testSeed1';       
testSeed1.strtype = 'r';
testSeed1.num     = [1 1];     
testSeed1.def     = @(val)nirs_get_defaults(...
    'readOnsets.addTestStimuli.testPType.testEP.NoFrequentSpikes.testSeed1', val{:}); 
testSeed1.help    = {'Enter a seed for the random number generator.'};

testSeed2         = cfg_entry; 
testSeed2.name    = 'Random seed';
testSeed2.tag     = 'testSeed2';       
testSeed2.strtype = 'r';
testSeed2.num     = [1 1];     
testSeed2.def     = @(val)nirs_get_defaults(...
    'readOnsets.addTestStimuli.testPType.testEP.FrequentSpikes.testSeed2', val{:}); 
testSeed2.help    = {'Enter a seed for the random number generator.'};
                        
testExpFastSpike         = cfg_entry; 
testExpFastSpike.name    = 'Time interval for frequent spikes';
testExpFastSpike.tag     = 'testExpFastSpike';       
testExpFastSpike.strtype = 'r';
testExpFastSpike.num     = [1 2];     
testExpFastSpike.def     = @(val)nirs_get_defaults(...
    'readOnsets.addTestStimuli.testPType.testEP.FrequentSpikes.testExpFastSpike', val{:}); 
testExpFastSpike.help    = {'Enter boundary min and max in seconds,'
                            'of the uniform distribution of time intervals,' 
                            'between frequent spikes.'}'; 

testExpSlowSpike1         = cfg_entry; 
testExpSlowSpike1.name    = 'Parameter for infrequent spikes interval';
testExpSlowSpike1.tag     = 'testExpSlowSpike1';       
testExpSlowSpike1.strtype = 'r';
testExpSlowSpike1.num     = [1 1];     
testExpSlowSpike1.def     = @(val)nirs_get_defaults(...
    'readOnsets.addTestStimuli.testPType.testEP.NoFrequentSpikes.testExpSlowSpike1', val{:}); 
testExpSlowSpike1.help    = {'Enter characteristic scale of infrequent spikes,'
                            'in seconds'}'; 

testExpSlowSpike2         = cfg_entry; 
testExpSlowSpike2.name    = 'Parameter for infrequent spikes interval';
testExpSlowSpike2.tag     = 'testExpSlowSpike2';       
testExpSlowSpike2.strtype = 'r';
testExpSlowSpike2.num     = [1 1];     
testExpSlowSpike2.def     = @(val)nirs_get_defaults(...
    'readOnsets.addTestStimuli.testPType.testEP.FrequentSpikes.testExpSlowSpike2', val{:}); 
testExpSlowSpike2.help    = {'Enter characteristic scale of infrequent spikes,'
                            'in seconds'}'; 

testAvgNumFastSpikes_perGroup         = cfg_entry; 
testAvgNumFastSpikes_perGroup.name    = 'Avg number of frequent spikes per group';
testAvgNumFastSpikes_perGroup.tag     = 'testAvgNumFastSpikes_perGroup';       
testAvgNumFastSpikes_perGroup.strtype = 'r';
testAvgNumFastSpikes_perGroup.num     = [1 1];     
testAvgNumFastSpikes_perGroup.def     = @(val)nirs_get_defaults(...
    'readOnsets.addTestStimuli.testPType.testEP.FrequentSpikes.testAvgNumFastSpikes_perGroup', val{:}); 
testAvgNumFastSpikes_perGroup.help    = {
            'Follows a Poisson process: Enter typical number.'
            'of fast spikes in a discharge'}'; 
        
testAvgNumSlowSpikes_perGroup         = cfg_entry; 
testAvgNumSlowSpikes_perGroup.name    = 'Average number of infrequent spikes per group';
testAvgNumSlowSpikes_perGroup.tag     = 'testAvgNumSlowSpikes_perGroup';       
testAvgNumSlowSpikes_perGroup.strtype = 'r';
testAvgNumSlowSpikes_perGroup.num     = [1 1];     
testAvgNumSlowSpikes_perGroup.def     = @(val)nirs_get_defaults(...
    'readOnsets.addTestStimuli.testPType.testEP.FrequentSpikes.testAvgNumSlowSpikes_perGroup', val{:}); 
testAvgNumSlowSpikes_perGroup.help    = {
            'Follows a Poisson process: Enter typical number.'
            'of infrequent spikes in a discharge'}'; 

testRescaleOn1         = cfg_menu;
testRescaleOn1.tag     = 'testRescaleOn1';
testRescaleOn1.name    = 'Rescale time to get specified number of spikes';
testRescaleOn1.help    = {'When rescaling, guaranteed to get the specified '
    'number of spikes. When not rescaling, the program will generate spikes '
    'until the specified number is reached, or the end of the dataset is reached, '
    'which could lead to a having fewer spikes than the number specified.'}';
testRescaleOn1.labels = {
               'True'
               'False'
                }';
testRescaleOn1.values = {1, 0};
testRescaleOn1.val = {1};

testRescaleOn2         = cfg_menu;
testRescaleOn2.tag     = 'testRescaleOn2';
testRescaleOn2.name    = 'Rescale time to get specified number of spikes';
testRescaleOn2.help    = {'When rescaling, guaranteed to get the specified '
    'number of spikes. When not rescaling, the program will generate spikes '
    'until the specified number is reached, or the end of the dataset is reached, '
    'which could lead to a having fewer spikes than the number specified.'}';
testRescaleOn2.labels = {
               'True'
               'False'
                }';
testRescaleOn2.values = {1, 0};
testRescaleOn2.val = {0};

NoFrequentSpikes           = cfg_branch; 
NoFrequentSpikes.name      = 'No frequent spikes';
NoFrequentSpikes.tag       = 'NoFrequentSpikes';
NoFrequentSpikes.val       = {testExpSlowSpike1 testRescaleOn1 testSeed1};
NoFrequentSpikes.help      = {'No frequent spikes, only one type of spikes'}; 
 
FrequentSpikes           = cfg_branch; 
FrequentSpikes.name      = 'Frequent spikes included';
FrequentSpikes.tag       = 'FrequentSpikes';
FrequentSpikes.val       = {testExpSlowSpike2 testExpFastSpike ...
     testAvgNumFastSpikes_perGroup testAvgNumSlowSpikes_perGroup testRescaleOn2 testSeed2};
FrequentSpikes.help      = {'Frequent spikes following a uniform distribution '
    'interspersed with slow spikes.'}'; 

%test Event paradigm
testEP           = cfg_choice; 
testEP.name      = 'Event Paradigm';
testEP.tag       = 'testEP';
testEP.values    = {NoFrequentSpikes FrequentSpikes};
testEP.val       = {FrequentSpikes};
%testEP.def     = @(val)nirs_get_defaults('readOnsets.addTestStimuli.testPType', val{:}); 
testEP.help      = {'Event Paradigm'
    'Choose whether to have only one type of spikes,'
    'all following an exponential distribution, '
    'or to also include fast spikes, which follow '
    'a uniform distribution.'
    'With two types of spikes, model works as follows,'
    'alternating between infrequent and frequent'
    'bunches of spikes (with numbers specified by separate Poisson distributions) '
    'putting the intervals one after the other. '
    'a number of infrequent spikes is '
    'obtained from a Poisson distribution; for each, a time interval to the next'
    'spike is generated from an exponential distribution '
    'This is then repeated for fast spikes, until the total desired number of' 
    'spikes is obtained. '}'; 

%Paradigm type
testPType      = cfg_choice;
testPType.tag  = 'testPType';
testPType.name = 'Paradigm: block or event';
testPType.values = {testBP testEP};
testPType.val = {testEP};
testPType.help = {'Choose test paradigm type: block (the number of stimuli '
    'will be evenly spread in time in the data file; or event (random stimuli '
    'generation, following a combination of uniform, exponential and Poisson distributions).'}';

%Test session number
testSessionNumber         = cfg_entry; 
testSessionNumber.name    = 'Test session number';
testSessionNumber.tag     = 'testSessionNumber';       
testSessionNumber.strtype = 'r';
testSessionNumber.num     = [1 1];     
testSessionNumber.def     = @(val)nirs_get_defaults(...
        'readOnsets.addTestStimuli.testSessionNumber', val{:}); 
testSessionNumber.help    = {'Enter NIRS data session number where stimuli are to be added.'}; 

%Test wavelength
testWavelength         = cfg_entry; 
testWavelength.name    = 'Test wavelength number(s)';
testWavelength.tag     = 'testWavelength';       
testWavelength.strtype = 'r';
testWavelength.num     = [1 Inf];     
testWavelength.def     = @(val)nirs_get_defaults(...
        'readOnsets.addTestStimuli.testWavelength', val{:}); 
testWavelength.help    = {'Enter wavelength number(s) of the data where the '
    'stimuli will be added. For example, enter 2 if OD at 690 nm is '
    'desired for the test and is the second wavelength.'}'; 

testAmplitude         = cfg_entry; 
testAmplitude.name    = 'Test response amplitude';
testAmplitude.tag     = 'testAmplitude';       
testAmplitude.strtype = 'r';
testAmplitude.num     = [1 Inf];     
testAmplitude.def     = @(val)nirs_get_defaults(...
            'readOnsets.addTestStimuli.testAmplitudeTarget.testAmplitude', val{:}); 
testAmplitude.help    = {'Enter desired response amplitude as a percentage. '
    'This will be understood as a percentage of the median of the raw intensity '
    'signal. If several channels are tested, a vector of percentages can '
    'be entered, that will be applied channelwise.'
    'effective SNR will be calculated for each channel and protocole.'}';

testSNR         = cfg_entry; 
testSNR.name    = 'Test SNR value';
testSNR.tag     = 'testSNR';       
testSNR.strtype = 'r';
testSNR.num     = [1 Inf];     
testSNR.def     = @(val)nirs_get_defaults(...
            'readOnsets.addTestStimuli.testAmplitudeTarget.testSNR', val{:}); 
testSNR.help    = {'Enter desired SNR in dB. '
    'The amplitude to be applied to each channel will be calculated '
    'Based on the formula SNR = 10 log10 (A * Ep/Eb),'
    'Where A is the amplitude, Eb is the power of the baseline (the sum'
    'of squares of amplitudes at each time point) and Ep is the power of'
    'the protocole after convolution with the HRFs.'
    'This option will result in a constant SNR value for each channel '
    'and each protocole, but the amplitude will vary.'
    'The calculated value for the amplitude will be stored.'}';

testAmplitudeTarget        = cfg_choice;
testAmplitudeTarget.name   = 'Target Amplitude or SNR';
testAmplitudeTarget.tag    = 'testAmplitudeTarget';
testAmplitudeTarget.values = {testAmplitude testSNR};
%Do not know how to specify a default value by a call using .def for a 
%cfg_choice object
testAmplitudeTarget.val    = {testSNR};
testAmplitudeTarget.help   = {'Choose whether target an amplitude level'
                        'Or a signal-to-noise ratio (SNR).'}';


testAmplitude2         = cfg_entry; 
testAmplitude2.name    = '2nd Volterra response amplitude';
testAmplitude2.tag     = 'testAmplitude2';       
testAmplitude2.strtype = 'r';
testAmplitude2.num     = [1 Inf];     
testAmplitude2.def     = @(val)nirs_get_defaults(...
            'readOnsets.addTestStimuli.testAmplitude2', val{:}); 
testAmplitude2.help    = {'Enter desired response amplitude as a percentage '
    'of the 1st Volterra kernel amplitude.'
    'Only used when 2nd Volterra is turned on.'}';

% ---------------------------------------------------------------------
% volt Model Interactions (Volterra)
% ---------------------------------------------------------------------
voltAddStim         = cfg_menu;
voltAddStim.tag     = 'voltAddStim';
voltAddStim.name    = 'Model Interactions (Volterra)';
voltAddStim.help    = {
                'Generalized convolution of inputs (U) with basis set (bf).'
                ''
                'For first order expansions the causes are simply convolved (e.g. stick functions) in U.u by the basis functions in bf to create a design matrix X.  For second order expansions new entries appear in ind, bf and name that correspond to the interaction among the orginal causes. The basis functions for these efects are two dimensional and are used to assemble the second order kernel. Second order effects are computed for only the first column of U.u.'
                'Interactions or response modulations can enter at two levels.  Firstly the stick function itself can be modulated by some parametric variate (this can be time or some trial-specific variate like reaction time) modeling the interaction between the trial and the variate or, secondly interactions among the trials themselves can be modeled using a Volterra series formulation that accommodates interactions over time (and therefore within and between trial types).'
                }';
voltAddStim.labels = {
               'Do not model Interactions'
               'Model Interactions'
                }';
voltAddStim.values = {1, 2};
voltAddStim.val = {1};

% Executable Branch
addTestStimuli      = cfg_exbranch;       
addTestStimuli.name = 'Add Stimuli with HRFs for testing';             
addTestStimuli.tag  = 'addTestStimuli'; 
addTestStimuli.val  = {NIRSmat DelPreviousData NewDirCopyNIRS testStimulusName testStimuliNumber ...
                testSessionNumber testWavelength testAmplitudeTarget ...
                voltAddStim testAmplitude2 keepAllChannels testChannels testDupChannels testPType};   
addTestStimuli.prog = @nirs_run_addTestStimuli;  
addTestStimuli.vout = @nirs_cfg_vout_addTestStimuli; 
addTestStimuli.help = {'Module to add stimuli convoluted with HRFs for'
            'simulation and testing purposes'}';

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_addTestStimuli(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration for preprocAnat: detect vitamins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Inputs:
% NIRSmat
% anatT1/image_in 

image_in         = cfg_files;
image_in.tag     = 'image_in';
image_in.name    = 'Image';
image_in.help    = {'Select the image to be processed. ROI will be selected in this image. Other stuff made with the image saved...'};
image_in.filter  = 'image';
image_in.ufilter = '.*';
image_in.num     = [1 1];

output_prefix_woVit         = cfg_entry;
output_prefix_woVit.name    = 'Prefix of the output anatomical image';
output_prefix_woVit.tag     = 'output_prefix_woVit';
output_prefix_woVit.strtype = 's';
output_prefix_woVit.num     = [1 Inf];
output_prefix_woVit.def  = @(val)nirs_get_defaults('preprocANAT.detectVitamins1.output_prefix_woVits', val{:});
output_prefix_woVit.help    = {'You can choose to give a particular prefix to the ',...
    'output image. This prefix will be added at the left of the name of the ',...
    'image. A default name will be given otherwise.'};


% Executable Branch
detectVitamins1      = cfg_exbranch;      
detectVitamins1.name = 'Coregistration with fiducials';            
detectVitamins1.tag  = 'detectVitamins1';
detectVitamins1.val  = {NIRSmat anatT1 output_prefix_woVit}; 
detectVitamins1.prog = @nirs_run_detectVitamins;  
detectVitamins1.vout = @nirs_cfg_vout_detectVitamins; 
detectVitamins1.help = {['This module detects fiducial markers (vitamin capsules)'...
    ' on anatomical image (T1), saves their positions in the NIRS.mat matrix, and'...
    ' creates a copy of the anatomical image where the markers are erased.']};

% Make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_detectVitamins(job)
vout = cfg_dep;                    
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration for MC segmentation: buildroi   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

image_in         = cfg_files;
image_in.tag     = 'image_in';
image_in.name    = 'Image';
image_in.help    = {'Select the image to be processed. ROI will be selected in this image. Other stuff made with the image saved...'};
image_in.filter  = 'image';
image_in.ufilter = '.*';
image_in.num     = [1 1];

% 
% crop_image      = cfg_menu;
% crop_image.tag  = 'crop_image';
% crop_image.name = 'Crop image';
% crop_image.labels = {'Yes','No'};
% crop_image.values = {1,0};
% crop_image.def  = @(val)nirs_get_defaults('preprocANAT.buildroi1.crop_image', val{:});
% crop_image.help = {'Help'};
    
%_______________________________________________________________________
output_prefix         = cfg_entry;
output_prefix.name    = 'Prefix of the output image';
output_prefix.tag     = 'output_prefix';
output_prefix.strtype = 's';
output_prefix.num     = [1 Inf];
output_prefix.def  = @(val)nirs_get_defaults('preprocANAT.buildroi1.output_prefix', val{:});
output_prefix.help    = {'You can choose to give a particular prefix to the ',...
    'output image. This prefix will be added at the left of the name of the ',...
    'image. A default name will be given otherwise.'};

%_______________________________________________________________________
buildroi1      = cfg_exbranch;
buildroi1.tag  = 'buildroi1';
buildroi1.name = 'Set vertices and build ROI';
buildroi1.val  = {NIRSmat keepAllChannels image_in output_prefix};
buildroi1.prog = @nirs_run_buildroi2;
buildroi1.help = {'Define region of interest.'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration for MC segmentation: MCsegment   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

image_in         = cfg_files;
image_in.tag     = 'image_in';
image_in.name    = 'Anatomical Images (Optional)';
image_in.help    = {'Select images to be processed. Same ROI will be selected ',...
     'for these images. Ancillary images will be saved... ',...
     'This module allows multi-subject processing, ',...
     'generating a segmented image for each subject. ',...
     'Note that a list of links to the segmented images will be ',...
     'available as a virtual output for further processing. ',...
     'It is essential that the order of the subjects be ',...
     'the same as in the readNIRS module.'};
image_in.filter  = 'image';
image_in.ufilter = '.*';
image_in.val{1}  = {''};
image_in.num     = [0 Inf];

%_______________________________________________________________________
output_autonaming      = cfg_menu;
output_autonaming.tag  = 'output_autonaming';
output_autonaming.name = 'Automatic output naming';
output_autonaming.labels = {'Yes','No'};
output_autonaming.values = {0,1};
output_autonaming.def = @(val)nirs_get_defaults('preprocANAT.MCsegment1.output_autonaming', val{:});
output_autonaming.help = {'Choose whether you want to choose the name of ',...
    'the output or not. If answer is ''Yes'', please enter name.'};
%_______________________________________________________________________
output_prefix      = cfg_entry;
output_prefix.tag  = 'output_prefix';
output_prefix.name = 'Prefix of the output image';
output_prefix.strtype = 's';
output_prefix.num     = [1 Inf];
output_prefix.def     = @(val)nirs_get_defaults('preprocANAT.MCsegment1.output_prefix', val{:});
output_prefix.help = {'You can choose to give a particular prefix to the ',...
    'output image. This prefix will be added at the left of the name of ',...
    'the image. A default name will be given otherwise.'};

%_______________________________________________________________________
sorting_method        = cfg_menu;
sorting_method.tag    = 'sorting_method';
sorting_method.name   = 'Sorting method';
sorting_method.labels = {'Median Filter','Opening','Gaussian Filter and Dilatation','Otsu','Median Filter and orthogonal Otsu'};
sorting_method.values = {0,1,2,3,4};
sorting_method.help   = {
    'NewSegment output images must be processed before pursuing segmentation : '
    '- (basic processing) Median Filter'
    '- Opening'
    '- Gaussian Filter and Dilatation'
    '- Otsu'
    'Explication de chacun des choix'
    }';

%_______________________________________________________________________
wtm      = cfg_branch;
wtm.tag  = 'wtm';
wtm.name = 'White matter';
sorting_method.def    = @(val)nirs_get_defaults('preprocANAT.MCsegment1.wtm', val{:});
wtm.val  = {sorting_method};
wtm.help = {'Options to process white matter image.'};
%_______________________________________________________________________
grm      = cfg_branch;
grm.tag  = 'grm';
grm.name = 'Grey matter';
sorting_method.def    = @(val)nirs_get_defaults('preprocANAT.MCsegment1.grm', val{:});
grm.val  = {sorting_method};
grm.help = {'Options to process grey matter image.'};
%_______________________________________________________________________
csf      = cfg_branch;
csf.tag  = 'csf';
csf.name = 'Cerebro-Spinal Fluid';
sorting_method.def    = @(val)nirs_get_defaults('preprocANAT.MCsegment1.csf', val{:});
csf.val  = {sorting_method};
csf.help = {'Options to process CSF image.'};
%_______________________________________________________________________
skl      = cfg_branch;
skl.tag  = 'skl';
skl.name = 'Skull';
sorting_method.def    = @(val)nirs_get_defaults('preprocANAT.MCsegment1.skl', val{:});
skl.val  = {sorting_method};
skl.help = {'Options to process skull image.'};
%_______________________________________________________________________
skn      = cfg_branch;
skn.tag  = 'skn';
skn.name = 'Skin';
sorting_method.def    = @(val)nirs_get_defaults('preprocANAT.MCsegment1.skn', val{:});
skn.val  = {sorting_method};
skn.help = {'Options to process skin image.'};

%_______________________________________________________________________
thresh_hs      = cfg_entry;
thresh_hs.tag  = 'thresh_hs';
thresh_hs.name = 'Threshold used to binarise head shadow in head shadow building';
thresh_hs.help = {'WHEN TO USE IT : if voxels are still zero after processing, you may increase this value keeping in mind, all boundaries will be blurred.'};
%_______________________________________________________________________
se_size_hs      = cfg_entry;
se_size_hs.tag  = 'se_size_hs';
se_size_hs.name = 'Size of the structural element in head shadow building';
se_size_hs.help = {
    'Size of the structural element (Mathematical morphology) used to clean the head shadow'
    'WHEN TO USE IT : (prefer change the threshold) nearly the same as the size of the threshold value.'
    }';
%_______________________________________________________________________
head_shadow      = cfg_branch;
head_shadow.tag  = 'head_shadow';
head_shadow.name = 'Build head shadowfor Monte Carlo Segmentation';

thresh_hs.def    = @(val)nirs_get_defaults('preprocANAT.head_shadow.thresh_hs', val{:});
se_size_hs.def   = @(val)nirs_get_defaults('preprocANAT.head_shadow.se_size_hs', val{:});
head_shadow.val  = {thresh_hs, se_size_hs};

head_shadow.help = {'The head shadow is obtained by summing all the ci-images (outputs from New Segment). It is used when sorting voxels inbetween one of the five layers. Once all the voxels have been sorted, some of them are zero whereas there are situated in the head volume, to solve this problem the head shadow is calculated and all the voxels belonging to the head are processed so that they are finally affected to one of the five layers'};

%_______________________________________________________________________
se_size_pi      = cfg_entry;
se_size_pi.tag  = 'se_size_pi';
se_size_pi.name = 'Size of the structural element in processing images';
se_size_pi.help = {
    'Size of the structural element (Mathematical morphology) used to clean images processed (ci-images ; outputs from New Segment)'
    'WHEN TO USE IT : .'
    }';
%_______________________________________________________________________
gaussfilt_size      = cfg_entry;
gaussfilt_size.tag  = 'gaussfilt_size';
gaussfilt_size.name = 'Size of the gaussian filter in processing images';
gaussfilt_size.help = {'...'
    'WHEN TO USE IT: if segmented image shows ''spatial instability''' 
    '(regions with a lot of different layers)'}';
%_______________________________________________________________________
gaussfilt_sdev      = cfg_entry;
gaussfilt_sdev.tag  = 'gaussfilt_sdev';
gaussfilt_sdev.name = 'Standard deviation of the gaussian filter in processing images';
gaussfilt_sdev.help = {
    '...'
    'WHEN TO USE IT : '};
%_______________________________________________________________________
process_image      = cfg_branch;
process_image.tag  = 'process_image';
process_image.name = 'Process image for Monte Carlo Segmentation';

se_size_pi.def     = @(val)nirs_get_defaults('preprocANAT.process_image.se_size_pi', val{:});
gaussfilt_size.def = @(val)nirs_get_defaults('preprocANAT.process_image.gaussfilt_size', val{:});
gaussfilt_sdev.def = @(val)nirs_get_defaults('preprocANAT.process_image.gaussfilt_sdev', val{:});
process_image.val  = {se_size_pi, gaussfilt_size, gaussfilt_sdev};

process_image.help = {'Gaussian Filter is only used in method Gaussian Filter and dilate.'};

%_______________________________________________________________________
thresh_as      = cfg_entry;
thresh_as.tag  = 'thresh_as';
thresh_as.name = 'Threshold for voxels with no conflict of belonging';
thresh_as.def         = @(val)nirs_get_defaults('preprocANAT.MCsegment1.thresh_as', val{:});
thresh_as.help = {'Once all the ci-images have been calculated by SPM, ',...
    'the algorithm sorts voxels with respect to the number of layers ',...
    'they belong to. For those belonging to only one layer, the ',...
    'algorithm verifies the random variable is near enough from 1 ',...
    '(bigger than the said threshold).'};

%_______________________________________________________________________
rebel_surrounding      = cfg_entry;
rebel_surrounding.tag  = 'rebel_surrounding';
rebel_surrounding.name = 'Surrounding size for rebel voxels';
rebel_surrounding.def = @(val)nirs_get_defaults('preprocANAT.MCsegment1.rebel_surrounding', val{:});
rebel_surrounding.help = {'Size of the surrounding used to processed and hence sort rebel voxels.'};

%_______________________________________________________________________
rebel_thresh_hs      = cfg_entry;
rebel_thresh_hs.tag  = 'rebel_thresh_hs';
rebel_thresh_hs.name = 'Threshold for head shadow in Monte Carlo segmentation';
rebel_thresh_hs.def   = @(val)nirs_get_defaults('preprocANAT.MCsegment1.rebel_thresh_hs', val{:});
rebel_thresh_hs.help = {'Threshold to discriminate voxels in head shadow image.'};

%_______________________________________________________________________
MCsegment1      = cfg_exbranch;
MCsegment1.tag  = 'MCsegment1';
MCsegment1.name = 'MC Segmentation';

MCsegment1.val  = {NIRSmat_optional DelPreviousData NewDirCopyNIRS image_in output_autonaming ...
    output_prefix skn skl csf grm wtm thresh_as head_shadow ...
    rebel_surrounding rebel_thresh_hs process_image};
MCsegment1.prog = @nirs_run_MCsegment;
MCsegment1.vout = @nirs_cfg_vout_MCsegment;
MCsegment1.help = {'Segmentation for Monte Carlo simulation and sorting ',...
    'of voxels of New Segment images. This module calls SPM''s New Segment ',...
    'unless it finds a c1 segmented file in the same directory as the input image.'};

%------------------------------------------------------------------------
function vout = nirs_cfg_vout_MCsegment(job)
% Determine what outputs will be present if this job is run. In this case,
% the structure of the inputs is fixed, and the output is always a single
% number. Note that input items may not be numbers, they can also be
% dependencies.
vout = cfg_dep;                     
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

% vout = cfg_dep;                        
% vout.sname      = 'Segmented Volume';    
% vout.src_output = substruct('.','Vsegmented'); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration for coregistration: coreg   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

segT1_4fit         = cfg_files; 
segT1_4fit.name    = 'Images segmented to fit P on.'; % The displayed name
segT1_4fit.tag     = 'segT1_4fit';       
segT1_4fit.filter  = 'image';    
segT1_4fit.num     = [0 Inf];
segT1_4fit.val{1}  = {''};
segT1_4fit.help    = {['(Optional) Choose the segmented image ',...
    '(0021_ or any other combination) to fit positions on scalp.',...
    'For multi-subject, the order of the images must correspond to ',...
    'the order of the subjects in the NIRS.mat matrix. ',...
    'If no input image is specified, the segmented image available in the ',...
    'NIRS matrix will be used (last one generated in MCsegment module).']};      

anatT1_template         = cfg_files; 
anatT1_template.name    = 'Anatomical template image'; 
anatT1_template.tag     = 'anatT1_template';       %file names
anatT1_template.filter  = 'image';  
anatT1_template.ufilter = '.*';
anatT1_template.def = @(val)nirs_get_defaults('coregNIRS.coreg1.anatT1_template', val{:});
anatT1_template.num     = [1 1];     % Number of inputs required 
anatT1_template.help    = {'Select anatomical template image for this subject.'}; 

fid_in_subject_MNI = cfg_menu;
fid_in_subject_MNI.tag    = 'fid_in_subject_MNI';
fid_in_subject_MNI.name   = 'Fiducials in subject MNI coordinates?';
fid_in_subject_MNI.labels = {'Yes','No'};
fid_in_subject_MNI.values = {1,0};
fid_in_subject_MNI.def    = @(val)nirs_get_defaults('coregNIRS.coreg1.fid_in_subject_MNI', val{:});
fid_in_subject_MNI.help   = {'Specify if coordinates of fiducials below are specified'
    'in the subject or patient MNI coordinates. The default option is NO: the default values'
    'are fiducial positions in normalized MNI coordinates of the SPM8 standard subject.'}';
        
%Atlas Normalized coordinates of fiducials
nasion_wMNI         = cfg_entry; %nasion_wMNI
nasion_wMNI.name    = 'Nasion';
nasion_wMNI.tag     = 'nasion_wMNI';       
nasion_wMNI.strtype = 'r';
nasion_wMNI.num     = [1 3];     
nasion_wMNI.def     = @(val)nirs_get_defaults('coregNIRS.coreg1.nasion_wMNI', val{:}); 
nasion_wMNI.help    = {'Coordinates of the nasion.'}; 

AL_wMNI         = cfg_entry; %AL_wMNI
AL_wMNI.name    = 'Auricular left';
AL_wMNI.tag     = 'AL_wMNI';       
AL_wMNI.strtype = 'r';
AL_wMNI.num     = [1 3];     
AL_wMNI.def     = @(val)nirs_get_defaults('coregNIRS.coreg1.AL_wMNI', val{:}); 
AL_wMNI.help    = {'Coordinates of Auricular Left.'}; 

AR_wMNI         = cfg_entry; %AR_wMNI
AR_wMNI.name    = 'Auricular right';
AR_wMNI.tag     = 'AR_wMNI';       
AR_wMNI.strtype = 'r';
AR_wMNI.num     = [1 3];     
AR_wMNI.def     = @(val)nirs_get_defaults('coregNIRS.coreg1.AR_wMNI', val{:}); 
AR_wMNI.help    = {'Coordinates of Auricular Right.'}; 

GenDataTopo = cfg_menu;
GenDataTopo.tag    = 'GenDataTopo';
GenDataTopo.name   = 'Generate topo data.';
GenDataTopo.labels = {'Yes','No'};
GenDataTopo.values = {1,0};
GenDataTopo.def    = @(val)nirs_get_defaults('coregNIRS.coreg1.GenDataTopo', val{:});
GenDataTopo.help   = {'Generate rend data (NIRS_SPM) for topographic '
            'reconstruction - stored in a separate file: TopoData.mat'}';


coreg1      = cfg_exbranch;       
coreg1.name = 'NIRScoreg';             
coreg1.tag  = 'coreg1'; 
coreg1.val  = {NIRSmat DelPreviousData NewDirCopyNIRS anatT1 segT1_4fit ...
    anatT1_template fid_in_subject_MNI nasion_wMNI AL_wMNI AR_wMNI GenDataTopo};    
coreg1.prog = @nirs_run_coreg;  
coreg1.vout = @nirs_cfg_vout_coreg; 
coreg1.help = {'Automatic coregistration.'};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_coreg(job)
vout = cfg_dep;                     
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration for coregistration: coreg MANUAL 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


NIRSmatSingle       = cfg_files; %Select NIRS.mat for this subject 
NIRSmatSingle.name    = 'Select NIRS.mat'; % The displayed name
NIRSmatSingle.tag     = 'NIRSmatSingle';       %file names
NIRSmatSingle.filter  = 'mat';
NIRSmatSingle.ufilter = '^NIRS.mat$';    
NIRSmatSingle.num     = [1 1];     % Number of inputs required 
NIRSmatSingle.help    = {'Select NIRS.mat for this subject.'}; 

CoregFromNIRS         = cfg_branch;
CoregFromNIRS.tag     = 'CoregFromNIRS';
CoregFromNIRS.name    = 'Manual coregistration using NIRS.mat';
CoregFromNIRS.val     = {NIRSmatSingle}; % tMCimg_config};
CoregFromNIRS.help    = {'Manual coregistration using previously created NIRS.mat'};

WanatT1         = cfg_files; %Select T1 for this subject 
WanatT1.name    = 'Select normalized anatomical image'; 
WanatT1.tag     = 'WanatT1';       %file names
WanatT1.filter  = 'image';  
WanatT1.ufilter = '.*';
WanatT1.num     = [1 1];     % Number of inputs required 
WanatT1.help    = {'Select normalized anatomical image for this subject.'}; 

NormParams         = cfg_files; %Select T1 for this subject 
NormParams.name    = 'Select normalized parameters'; 
NormParams.tag     = 'NormParams';       %file names
%NormParams.filter  = '';  
NormParams.ufilter = '_sn.*';
NormParams.num     = [1 1];     % Number of inputs required 
NormParams.help    = {'Select normalization parameters for this subject.'}; 


Coreg_standalone         = cfg_branch;
Coreg_standalone.tag     = 'Coreg_standalone';
Coreg_standalone.name    = 'Standalone coregistration';
Coreg_standalone.val     = {WanatT1 NormParams}; % tMCimg_config};
Coreg_standalone.help    = {'Select normalized anatomical image for optode positioning for this subject'};

Coreg_choice        = cfg_choice;
Coreg_choice.name   = 'Manual coregistration choice';
Coreg_choice.tag    = 'Coreg_choice';
Coreg_choice.values = {CoregFromNIRS,Coreg_standalone};
Coreg_choice.val    = {CoregFromNIRS};
Coreg_choice.help   = {'Choose whether to run from NIRS.mat or as standalone'};

Vsegmented         = cfg_files; %Select anatomical image for this subject 
Vsegmented.name    = 'Select anatomical image'; % The displayed name
Vsegmented.tag     = 'Vsegmented';       %file names
Vsegmented.filter  = 'image';  
Vsegmented.ufilter = '.*';
Vsegmented.num     = [1 1];     % Number of inputs required 
Vsegmented.help    = {'Select native (not normalized) anatomical image for this subject.'};

% Executable Branch
coreg_manual1      = cfg_exbranch;      
coreg_manual1.name = 'Manual NIRScoreg';            
coreg_manual1.tag  = 'coreg_manual1';
coreg_manual1.val  = {Coreg_choice Vsegmented }; 
coreg_manual1.prog = @nirs_run_coreg_manual;  
coreg_manual1.vout = @nirs_cfg_vout_coreg_manual; 
coreg_manual1.help = {'Manual coregistration.'};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_coreg_manual(job)
vout = cfg_dep;                    
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration for coregistration: view 3D  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

on_cortex        = cfg_menu;
on_cortex.tag    = 'on_cortex';
on_cortex.name   = 'View optodes on cortex.';
on_cortex.labels = {'Yes','No'};
on_cortex.values = {1,0};
on_cortex.def    = @(val)nirs_get_defaults('coregNIRS.view3d1.on_cortex', val{:});
on_cortex.help   = {'If answer is No, default view consists on displaying optodes on the scalp.'};

save_figure        = cfg_menu;
save_figure.tag    = 'save_figure';
save_figure.name   = 'Save 3D view figure.';
save_figure.labels = {'Yes','No'};
save_figure.values = {1,0};
save_figure.def    = @(val)nirs_get_defaults('coregNIRS.view3d1.save_figure', val{:});
save_figure.help   = {'Yes: 3D view will be saved as Matlab figure; No: figure will remain open.'};

segT1_4fit         = cfg_files; 
segT1_4fit.name    = 'Images segmented viewed in 3D.'; % The displayed name
segT1_4fit.tag     = 'segT1_4fit';       
segT1_4fit.filter  = 'image';    
segT1_4fit.num     = [0 Inf];   
segT1_4fit.val{1}  = {''};
segT1_4fit.help    = {['Optodes will be represented on the cortex or the ',...
    'scalp depending on the previous answer. If you answered Yes, please ',...
    'choose ''c3_'' image, otherwise choose the segmented image ',...
    '(0021_ or any other combination).',...
    'For multi-subject, the order of the images must correspond to ',...
    'the order of the subjects in the NIRS.mat matrix. ',...
    'If no input image is specified, the segmented image available in the ',...
    'NIRS matrix will be used (last one generated in MCsegment module).']};        

% Executable Branch
view3d1      = cfg_exbranch;       % This is the branch that has information about how to run this module
view3d1.name = '3D view to check coregistration';             % The display name
view3d1.tag  = 'view3d1'; %Very important: tag is used when calling for execution
view3d1.val  = {NIRSmat on_cortex segT1_4fit save_figure};
view3d1.prog = @nirs_run_view3d;  % A function handle that will be called with the harvested job to run the computation
%view3d1.vout = @nirs_cfg_vout_view3d; % A function handle that will be called with the harvested job to determine virtual outputs
view3d1.help = {'Display 3D view of the skin with optodes so that user is able to check coregistration.'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Module 4: Utilities for NIRS data preprocessing: heart rate, filters, etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4.0 Paces... Heart and Mayer -- Mayer not done yet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

threshold_stdev         = cfg_entry;
threshold_stdev.name    = 'Threshold on standard deviation';
threshold_stdev.tag     = 'threshold_stdev';
threshold_stdev.strtype = 'r';
threshold_stdev.num     = [1 1];
threshold_stdev.def     = @(val)nirs_get_defaults('preprocessNIRS.remove_chn_stdev.threshold_stdev', val{:});
threshold_stdev.help    = {'Enter cutoff as a percentage (relative to median of channel) ',...
    'of the median of the standard deviation calculated in rolling windows'.'};

window_stdev         = cfg_entry;
window_stdev.name    = 'Rolling window length';
window_stdev.tag     = 'window_stdev';
window_stdev.strtype = 'r';
window_stdev.num     = [1 1];
window_stdev.def     = @(val)nirs_get_defaults('preprocessNIRS.remove_chn_stdev.window_stdev', val{:});
window_stdev.help    = {'Enter the length in seconds for the rolling windows.'};

% Executable Branch
remove_chn_stdev      = cfg_exbranch;       
remove_chn_stdev.name = 'Remove noisy channels (stdev)';             
remove_chn_stdev.tag  = 'remove_chn_stdev'; 
remove_chn_stdev.val  = {NIRSmat DelPreviousData NewDirCopyNIRS threshold_stdev window_stdev};   
remove_chn_stdev.prog = @nirs_run_remove_chn_stdev;  
remove_chn_stdev.vout = @nirs_cfg_vout_remove_chn_stdev;
remove_chn_stdev.help = {['Preprocessing step: remove noisy channels ',...
    'on the median of a channelwise rolling standard deviation measure.']};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_remove_chn_stdev(job)
vout = cfg_dep;                    
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4.1 Paces... Heart and Mayer -- Mayer not done yet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Remove channels without a clear heart beat
remove_no_heartbeat      = cfg_menu;
remove_no_heartbeat.tag  = 'remove_no_heartbeat';
remove_no_heartbeat.name = 'Remove Channels without clear heart beat';
remove_no_heartbeat.labels = {'True','False'};
remove_no_heartbeat.values = {1,0};
remove_no_heartbeat.def  = @(val)nirs_get_defaults('preprocessNIRS.criugm_paces1.remove_no_heartbeat', val{:});
remove_no_heartbeat.help = {['Remove channels without clear heart beat; ',...
        'Detection carried out only on wavelength most sensitive to heart beat. ',...
        'If no heart beat found at that wavelength, the other wavelengths are removed too.'] };

%parameters for heart_resting
win_type = cfg_menu;
win_type.tag  = 'win_type';
win_type.name = 'Type of window for probes';
win_type.labels = {'Hanning','Hamming','Rect'};
win_type.values = {0,1,2};
win_type.def  = @(val)nirs_get_defaults(...
    'preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_resting.STFT_param.win_type', val{:});
win_type.help = {'Only Hanning for now.'};

win_width         = cfg_entry;
win_width.name    = 'Window width';
win_width.tag     = 'win_width';       
win_width.strtype = 'r';
win_width.num     = [1 1];
win_width.def  = @(val)nirs_get_defaults(...
    'preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_resting.STFT_param.win_width', val{:});
win_width.help    = {'Window width in SECONDS'};

Nprobe         = cfg_entry;
Nprobe.name    = 'Number of probes';
Nprobe.tag     = 'Nprobe';       
Nprobe.strtype = 'r';
Nprobe.num     = [1 1];
Nprobe.def  = @(val)nirs_get_defaults(...
    'preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_resting.STFT_param.Nprobe', val{:});
Nprobe.help    = {'Number of probes taken along the signal (power of 2)'};

fft_size         = cfg_entry;
fft_size.name    = 'FFT size';
fft_size.tag     = 'fft_size';       
fft_size.strtype = 'r';
fft_size.num     = [1 1];
fft_size.def  = @(val)nirs_get_defaults(...
    'preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_resting.STFT_param.fft_size', val{:});
fft_size.help    = {'FFT size for each probe'};

STFT_param         = cfg_branch;
STFT_param.tag     = 'STFT_param';
STFT_param.name    = 'Parameters for heart pace search';
STFT_param.val     = {win_type win_width Nprobe fft_size};
STFT_param.help    = {'Short Term Fourier Transform Parameters'};

%Detection wavelengths
% detect_wavelength         = cfg_entry; 
% detect_wavelength.name    = 'Detection wavelength number';
% detect_wavelength.tag     = 'detect_wavelength';       
% detect_wavelength.strtype = 'r';
% detect_wavelength.num     = [1 Inf];     
% detect_wavelength.def     = @(val)nirs_get_defaults(...
%     'preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_resting.detect_wavelength', val{:}); 
% detect_wavelength.help    = {['Enter wavelength number(s) for detection of ',...
%     'the heart rate. If no heart rate is detected, and the remove channels ',...
%     'option is selected, channels for all wavelengths at this location will be removed. ',...
%     'For example, enter 1 if OD at 830 nm is the first wavelength; enter an array, ',...
%     'such as [1 2] if detection of heart rate at the first two wavelengths is required.']}; 
% 
MinHeartRate         = cfg_entry; 
MinHeartRate.name    = 'Minimum Heart Rate for Detection';
MinHeartRate.tag     = 'MinHeartRate';       
MinHeartRate.strtype = 'r';
MinHeartRate.num     = [1 1];     
MinHeartRate.def     = @(val)nirs_get_defaults(...
    'preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_resting.MinHeartRate', val{:}); 
MinHeartRate.help    = {['Enter minimum heart rate allowed for final detection in beats per minute.']};

MaxHeartRate         = cfg_entry; 
MaxHeartRate.name    = 'Maximum Heart Rate for Detection';
MaxHeartRate.tag     = 'MaxHeartRate';       
MaxHeartRate.strtype = 'r';
MaxHeartRate.num     = [1 1];     
MaxHeartRate.def     = @(val)nirs_get_defaults(...
    'preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_resting.MaxHeartRate', val{:}); 
MaxHeartRate.help    = {['Enter maximum heart rate allowed for final detection in beats per minute.']};

InternalMinHeartRate         = cfg_entry; 
InternalMinHeartRate.name    = 'Internal Minimum Heart Rate for Detection';
InternalMinHeartRate.tag     = 'InternalMinHeartRate';       
InternalMinHeartRate.strtype = 'r';
InternalMinHeartRate.num     = [1 1];     
InternalMinHeartRate.def     = @(val)nirs_get_defaults(...
    'preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_resting.InternalMinHeartRate', val{:}); 
InternalMinHeartRate.help    = {['Enter minimum heart rate allowed for detection in beats per minute, ',...
                        '"internal", i.e. for check over small data windows for FFT.']};

InternalMaxHeartRate         = cfg_entry; 
InternalMaxHeartRate.name    = 'Internal Maximum Heart Rate for Detection';
InternalMaxHeartRate.tag     = 'InternalMaxHeartRate';       
InternalMaxHeartRate.strtype = 'r';
InternalMaxHeartRate.num     = [1 1];     
InternalMaxHeartRate.def     = @(val)nirs_get_defaults(...
    'preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_resting.InternalMaxHeartRate', val{:}); 
InternalMaxHeartRate.help    = {['Enter maximum heart rate allowed for detection in beats per minute, ',...
                        '"internal", i.e. for check over small data windows for FFT.']};

MaxHeartStdev         = cfg_entry; 
MaxHeartStdev.name    = 'Maximum Heart Rate Standard Deviation';
MaxHeartStdev.tag     = 'MaxHeartStdev';       
MaxHeartStdev.strtype = 'r';
MaxHeartStdev.num     = [1 1];     
MaxHeartStdev.def     = @(val)nirs_get_defaults(...
    'preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_resting.MaxHeartStdev', val{:}); 
MaxHeartStdev.help    = {'Enter maximum heart rate standard deviation allowed in beats per minute.'}';

%copy all parameters with "2" for heart_exercise
win_type2 = cfg_menu;
win_type2.tag  = 'win_type2';
win_type2.name = 'Type of window for probes';
win_type2.labels = {'Hanning','Hamming','Rect'};
win_type2.values = {0,1,2};
win_type2.def  = @(val)nirs_get_defaults(...
    'nirs10.preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_exercise.STFT_param2.win_type2', val{:});
win_type2.help = {'Only Hanning for now.'};

win_width2         = cfg_entry;
win_width2.name    = 'Window width';
win_width2.tag     = 'win_width2';       
win_width2.strtype = 'r';
win_width2.num     = [1 1];
win_width2.def  = @(val)nirs_get_defaults(...
    'nirs10.preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_exercise.STFT_param2.win_width2', val{:});
win_width2.help    = {'Window width in SECONDS'};

Nprobe2         = cfg_entry;
Nprobe2.name    = 'Number of probes';
Nprobe2.tag     = 'Nprobe2';       
Nprobe2.strtype = 'r';
Nprobe2.num     = [1 1];
Nprobe2.def  = @(val)nirs_get_defaults(...
    'nirs10.preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_exercise.STFT_param2.Nprobe2', val{:});
Nprobe2.help    = {'Number of probes taken along the signal (power of 2)'};

fft_size2         = cfg_entry;
fft_size2.name    = 'FFT size';
fft_size2.tag     = 'fft_size2';       
fft_size2.strtype = 'r';
fft_size2.num     = [1 1];
fft_size2.def  = @(val)nirs_get_defaults(...
    'nirs10.preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_exercise.STFT_param2.fft_size2', val{:});
fft_size2.help    = {'FFT size for each probe'};

STFT_param2         = cfg_branch;
STFT_param2.tag     = 'STFT_param2';
STFT_param2.name    = 'Parameters for heart pace search';
STFT_param2.val     = {win_type2 win_width2 Nprobe2 fft_size2};
STFT_param2.help    = {'Short Term Fourier Transform Parameters'};

%Detection wavelengths
% detect_wavelength2         = cfg_entry; 
% detect_wavelength2.name    = 'Detection wavelength number';
% detect_wavelength2.tag     = 'detect_wavelength2';       
% detect_wavelength2.strtype = 'r';
% detect_wavelength2.num     = [1 Inf];     
% detect_wavelength2.def     = @(val)nirs_get_defaults(...
%     'nirs10.preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_exercise.detect_wavelength2', val{:}); 
% detect_wavelength2.help    = {['Enter wavelength number(s) for detection of ',...
%     'the heart rate. If no heart rate is detected, and the remove channels ',...
%     'option is selected, channels for all wavelengths at this location will be removed. ',...
%     'For example, enter 1 if OD at 830 nm is the first wavelength; enter an array, ',...
%     'such as [1 2] if detection of heart rate at the first two wavelengths is required.']}; 

% MinHeartRate2         = cfg_entry; 
% MinHeartRate2.name    = 'Minimum Heart Rate for Detection';
% MinHeartRate2.tag     = 'MinHeartRate2';       
% MinHeartRate2.strtype = 'r';
% MinHeartRate2.num     = [1 1];     
% MinHeartRate2.def     = @(val)nirs_get_defaults(...
%     'nirs10.preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_exercise.MinHeartRate2', val{:}); 
% MinHeartRate2.help    = {'Enter minimum heart rate allowed for final detection in Hz.'}';
% 
% MaxHeartRate2         = cfg_entry; 
% MaxHeartRate2.name    = 'Maximum Heart Rate for Detection';
% MaxHeartRate2.tag     = 'MaxHeartRate2';       
% MaxHeartRate2.strtype = 'r';
% MaxHeartRate2.num     = [1 1];     
% MaxHeartRate2.def     = @(val)nirs_get_defaults(...
%     'nirs10.preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_exercise.MaxHeartRate2', val{:}); 
% MaxHeartRate2.help    = {'Enter maximum heart rate allowed for final detection in Hz.'}';
% 
% InternalMinHeartRate2         = cfg_entry; 
% InternalMinHeartRate2.name    = 'Internal Minimum Heart Rate for Detection';
% InternalMinHeartRate2.tag     = 'InternalMinHeartRate2';       
% InternalMinHeartRate2.strtype = 'r';
% InternalMinHeartRate2.num     = [1 1];     
% InternalMinHeartRate2.def     = @(val)nirs_get_defaults(...
%     'nirs10.preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_exercise.InternalMinHeartRate2', val{:}); 
% InternalMinHeartRate2.help    = {'Enter minimum heart rate allowed for detection in Hz, '
%                         '"internal", i.e. for check over small data windows for FFT.'}';
% 
% InternalMaxHeartRate2         = cfg_entry; 
% InternalMaxHeartRate2.name    = 'Internal Maximum Heart Rate for Detection';
% InternalMaxHeartRate2.tag     = 'InternalMaxHeartRate2';       
% InternalMaxHeartRate2.strtype = 'r';
% InternalMaxHeartRate2.num     = [1 1];     
% InternalMaxHeartRate2.def     = @(val)nirs_get_defaults(...
%     'nirs10.preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_exercise.InternalMaxHeartRate2', val{:}); 
% InternalMaxHeartRate2.help    = {'Enter maximum heart rate allowed for detection in Hz, '
%                         '"internal", i.e. for check over small data windows for FFT.'}';
% 
% MaxHeartStdev2         = cfg_entry; 
% MaxHeartStdev2.name    = 'Maximum Heart Rate Standard Deviation';
% MaxHeartStdev2.tag     = 'MaxHeartStdev2';       
% MaxHeartStdev2.strtype = 'r';
% MaxHeartStdev2.num     = [1 1];     
% MaxHeartStdev2.def     = @(val)nirs_get_defaults(...
%     'nirs10.preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_exercise.MaxHeartStdev2', val{:}); 
% MaxHeartStdev2.help    = {'Enter maximum heart rate standard deviation allowed in beats per second'}';

heart_resting         = cfg_branch;
heart_resting.tag     = 'heart_resting';
heart_resting.name    = 'Resting state parameters';
heart_resting.val     = {STFT_param MinHeartRate MaxHeartRate ...
    InternalMinHeartRate InternalMaxHeartRate MaxHeartStdev}; 
heart_resting.help    = {'Choose parameters for resting state heart rate detection'};

heart_exercise         = cfg_branch;
heart_exercise.tag     = 'heart_exercise';
heart_exercise.name    = 'Parameters during exercise';
heart_exercise.val     = {STFT_param}; %detect_wavelength MinHeartRate MaxHeartRate ...
%     InternalMinHeartRate InternalMaxHeartRate MaxHeartStdev}; 
heart_exercise.help    = {'Choose parameters for heart rate detection'
                        'during aerobic exercise (such as VO2max test).'}';

heart_rate_cfg           = cfg_choice;
heart_rate_cfg.name      = 'Heart rate configuration';
heart_rate_cfg.tag       = 'heart_rate_cfg';
heart_rate_cfg.values    = {heart_resting heart_exercise};
%heart_rate_cfg.def     =
%@(val)nirs_get_defaults('preprocessNIRS.criugm_paces1', val{:}); 
heart_rate_cfg.val       = {heart_resting}; 
heart_rate_cfg.help      = {'Choose configuration of parameters.'
                            'Resting-state or VO2max (aerobic exercise)'}'; 

% Executable Branch
criugm_paces1      = cfg_exbranch;       
criugm_paces1.name = 'Heart rate utility';             
criugm_paces1.tag  = 'criugm_paces1'; 
criugm_paces1.val  = {NIRSmat DelPreviousData NewDirCopyNIRS heart_rate_cfg remove_no_heartbeat};   
criugm_paces1.prog = @nirs_run_criugm_paces;  
criugm_paces1.vout = @nirs_cfg_vout_criugm_paces;
criugm_paces1.help = {['Preprocessing step: Extract heart rate and, if desired, ',...
    'remove channels without a clear detectable heart rate.']}';

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_criugm_paces(job)
vout = cfg_dep;                     
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Preprocess NIRS: mark negative points as bad and interpolate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sum_neg_threshold         = cfg_entry;
sum_neg_threshold.name    = 'Threshold on number of channels';
sum_neg_threshold.tag     = 'sum_neg_threshold';
sum_neg_threshold.strtype = 'r';
sum_neg_threshold.num     = [1 1];
sum_neg_threshold.def     = @(val)nirs_get_defaults('preprocessNIRS.mark_negative.sum_neg_threshold', val{:});
sum_neg_threshold.help    = {['Enter the value of the threshold as ',...
    'a percentage of the total number of channels, so that ',...
    'data points which are marked bad for at least that many channels ',...
    'will be marked bad for all channels.']};

% Executable Branch
mark_negative      = cfg_exbranch;       
mark_negative.name = 'Mark negative and interpolate';             
mark_negative.tag  = 'mark_negative'; 
mark_negative.val  = {NIRSmat DelPreviousData NewDirCopyNIRS sum_neg_threshold};   
mark_negative.prog = @nirs_run_mark_negative;  
mark_negative.vout = @nirs_cfg_vout_mark_negative;
mark_negative.help = {['Preprocessing step: mark negative data points as ',...
    'bad data points, and interpolate from nearby values.']};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_mark_negative(job)
vout = cfg_dep;                    
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Preprocess NIRS: mark movement jumps and shifts as bad and normalize 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mvt_ch_thresh         = cfg_entry;
mvt_ch_thresh.name    = 'Movement time cutoff';
mvt_ch_thresh.tag     = 'mvt_ch_thresh';
mvt_ch_thresh.strtype = 'r';
mvt_ch_thresh.num     = [1 1];
mvt_ch_thresh.def     = @(val)nirs_get_defaults('preprocessNIRS.mark_movement.mvt_ch_thresh', val{:});
mvt_ch_thresh.help    = {'Enter the maximal percentage of time allowed '
    'for a good channel - this will be tested in the first session'
    'channels that failed will be excluded for this and following sessions.'}';

mvt_window_length         = cfg_entry;
mvt_window_length.name    = 'Window size to detect movement';
mvt_window_length.tag     = 'mvt_window_length';
mvt_window_length.strtype = 'r';
mvt_window_length.num     = [1 1];
mvt_window_length.def     = @(val)nirs_get_defaults('preprocessNIRS.mark_movement.mvt_window_length', val{:});
mvt_window_length.help    = {['Enter the length of the window in seconds ',...
    'over which to detect movement.']};

mvt_cutoff         = cfg_entry;
mvt_cutoff.name    = 'Movement cutoff';
mvt_cutoff.tag     = 'mvt_cutoff';
mvt_cutoff.strtype = 'r';
mvt_cutoff.num     = [1 1];
mvt_cutoff.def     = @(val)nirs_get_defaults('preprocessNIRS.mark_movement.mvt_cutoff', val{:});
mvt_cutoff.help    = {'Enter the maximal change allowed '
    'as a percentage of the median intensity signal over the window '
    'length, above which data points over the window length are marked as bad.'}';

sum_mvt_threshold         = cfg_entry;
sum_mvt_threshold.name    = 'Threshold on number of channels';
sum_mvt_threshold.tag     = 'sum_mvt_threshold';
sum_mvt_threshold.strtype = 'r';
sum_mvt_threshold.num     = [1 1];
sum_mvt_threshold.def     = @(val)nirs_get_defaults('preprocessNIRS.mark_movement.sum_mvt_threshold', val{:});
sum_mvt_threshold.help    = {'Enter the value of the threshold as '
    'a percentage of the total number of channels, so that '
    'data points which are marked bad for at least that many channels '
    'will be marked bad for all channels.'}';

min_session_duration         = cfg_entry;
min_session_duration.name    = 'Minimum length of subsessions';
min_session_duration.tag     = 'min_session_duration';
min_session_duration.strtype = 'r';
min_session_duration.num     = [1 1];
min_session_duration.def     = @(val)nirs_get_defaults('preprocessNIRS.mark_movement.min_session_duration', val{:});
min_session_duration.help    = {'Enter the minimum length (for example 60 seconds).'
                        'of subsessions, i.e. the minimum intervals between '
                        'markers of movement.'}';

% Executable Branch
mark_movement      = cfg_exbranch;       
mark_movement.name = 'Mark movement';             
mark_movement.tag  = 'mark_movement'; 
mark_movement.val  = {NIRSmat DelPreviousData NewDirCopyNIRS mvt_ch_thresh...
    mvt_window_length mvt_cutoff sum_mvt_threshold min_session_duration};   
mark_movement.prog = @nirs_run_mark_movement;  
mark_movement.vout = @nirs_cfg_vout_mark_movement;
mark_movement.help = {['Preprocessing step: mark movement jumps ',...
    'and shifts, segment into intervals.']};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_mark_movement(job)
vout = cfg_dep;                    
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.4. Normalize to baseline  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Choose to normalize Optical Densities to initial value or to median    
Normalize_OD = cfg_menu;
Normalize_OD.tag  = 'Normalize_OD';
Normalize_OD.name = 'Normalization method';
Normalize_OD.labels = {'Median','Initial Value','Mean'};
Normalize_OD.values = {0,1,2};
Normalize_OD.def  = @(val)nirs_get_defaults('preprocessNIRS.normalize_baseline.Normalize_OD', val{:});
Normalize_OD.help = {['Choose normalization of Optical Densities',...
        'Median (preferred), Initial value, or Mean.']};
    
add_or_mult      = cfg_menu;
add_or_mult.tag  = 'add_or_mult';
add_or_mult.name = 'Additive or multiplicative';
add_or_mult.labels = {'Additive', 'Multiplicative'};
add_or_mult.values = {1,0};
add_or_mult.def  = @(val)nirs_get_defaults('preprocessNIRS.normalize_baseline.add_or_mult', val{:});
add_or_mult.help = {'Select whether using additive (on concentrations for example)'
    'or multiplicative (on optical intensities) normalization.' }';

baseline_duration         = cfg_entry; 
baseline_duration.name    = 'Baseline duration'; % The displayed name
baseline_duration.tag     = 'baseline_duration';       %file names
baseline_duration.strtype = 'r';  
baseline_duration.num     = [1 1];     % Number of inputs required 
%baseline_duration.val{1}  = 100;
baseline_duration.def     = @(val)nirs_get_defaults('preprocessNIRS.normalize_baseline.baseline_duration', val{:});
baseline_duration.help    = {'Enter the baseline duration in seconds to use '
                    'prior to stimuli - applies only '
                    'to the normalization type by stimuli below.)'}'; 

normalization_type      = cfg_menu;
normalization_type.tag  = 'normalization_type';
normalization_type.name = 'Normalization type';
normalization_type.labels = {'Global', 'By bad point segments', 'By stimuli'};
normalization_type.values = {1,2,3};
normalization_type.def  = @(val)nirs_get_defaults('preprocessNIRS.normalize_baseline.normalization_type', val{:});
normalization_type.help = {'Normalization type: global, by bad point segments, before stimuli.'
    'When normalizing by stimuli, the code finds each stimulus instance '
    'and uses the time window specified prior to the stimulus for normalization.'}'; 

Analyzer_sf         = cfg_entry; 
Analyzer_sf.name    = 'Scaling factor'; % The displayed name
Analyzer_sf.tag     = 'Analyzer_sf';       %file names
Analyzer_sf.strtype = 'r';  
Analyzer_sf.num     = [1 1];     % Number of inputs required 
%Analyzer_sf.val{1}  = 100;
Analyzer_sf.def     = @(val)nirs_get_defaults('preprocessNIRS.normalize_baseline.Analyzer_sf', val{:});
Analyzer_sf.help    = {'Apply a scaling factor on the amplitude of '
                    'all channels (for easier visualization with Analyzer.)'}'; 
                
% Executable Branch
normalize_baseline      = cfg_exbranch;       
normalize_baseline.name = 'Normalize Baseline';             
normalize_baseline.tag  = 'normalize_baseline'; 
normalize_baseline.val  = {NIRSmat DelPreviousData NewDirCopyNIRS Normalize_OD add_or_mult ...
        baseline_duration normalization_type Analyzer_sf};   
normalize_baseline.prog = @nirs_run_normalize_baseline;  
normalize_baseline.vout = @nirs_cfg_vout_normalize_baseline;
normalize_baseline.help = {'Normalize to baseline'}';

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_normalize_baseline(job)
vout = cfg_dep;                    
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration for converting Optical Densities to HbO/HbR  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PVF         = cfg_entry;
PVF.name    = 'Partial Volume Factors';
PVF.tag     = 'PVF';
PVF.strtype = 'r';
PVF.num     = [1 2];
PVF.def  = @(val)nirs_get_defaults('preprocessNIRS.ODtoHbOHbR.PVF', val{:});
PVF.help    = {'Enter the partial volume factor values for each wavelength ',...
    'as a vector: [PVF(lambda_1) ... PVF(lambda_n)].'};


% ---------------------------------------------------------------------
% lpf Low-pass filter
% ---------------------------------------------------------------------
fwhm2      = cfg_entry;
fwhm2.tag  = 'fwhm2';
fwhm2.name = 'FWHM in seconds';
fwhm2.val = {1.5};
fwhm2.strtype = 'r';  
fwhm2.num     = [1 1]; 
fwhm2.help    = {'FWHM in seconds.'}; 

downsamplingFactor      = cfg_entry;
downsamplingFactor.tag  = 'downsamplingFactor';
downsamplingFactor.name = 'Downsampling Factor';
downsamplingFactor.val = {1}; %10
downsamplingFactor.strtype = 'r';  
downsamplingFactor.num     = [1 1]; 
downsamplingFactor.help    = {'Specify downsampling factor.'}; 

downsampleWhen         = cfg_menu;
downsampleWhen.tag     = 'downsampleWhen';
downsampleWhen.name    = 'Apply downsampling on';
downsampleWhen.help    = {'Choose at what step to apply low-pass filtering and downsampling.'};
downsampleWhen.labels = {
                 'Raw Optical Densities'
                 'Log of Optical Densities'
                 'Concentrations'}';
downsampleWhen.values = {1 2 3};
downsampleWhen.val    = {1};

lpf_gauss2         = cfg_branch;
lpf_gauss2.tag     = 'lpf_gauss2';
lpf_gauss2.name    = 'Gaussian Filter';
lpf_gauss2.val     = {fwhm2 downsamplingFactor downsampleWhen}; 
lpf_gauss2.help    = {'Specify properties of Gaussian filter'};

lpf_none         = cfg_branch;
lpf_none.tag     = 'lpf_none';
lpf_none.name    = 'No low pass filter';
lpf_none.help    = {'No low pass filter.'};

nirs_lpf2           = cfg_choice;
nirs_lpf2.name      = 'Low-pass filter';
nirs_lpf2.tag       = 'nirs_lpf2';
nirs_lpf2.values    = {lpf_none
                      lpf_gauss2}; 
nirs_lpf2.val       = {lpf_none}; 
nirs_lpf2.help      = {'Choose low-pass filter.'}; 


% Executable Branch
ODtoHbOHbR      = cfg_exbranch;       
ODtoHbOHbR.name = 'Convert OD to HbO/HbR ';             
ODtoHbOHbR.tag  = 'ODtoHbOHbR'; 
ODtoHbOHbR.val  = {NIRSmat DelPreviousData NewDirCopyNIRS PVF nirs_lpf2}; 
ODtoHbOHbR.prog = @nirs_run_ODtoHbOHbR;  
ODtoHbOHbR.vout = @nirs_cfg_vout_ODtoHbOHbR; 
ODtoHbOHbR.help = {'Convert OD to HbO/HbR.'}';

function vout = nirs_cfg_vout_ODtoHbOHbR(job)
vout = cfg_dep;                     
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.6 High pass and low pass filters (optional)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------------------------------------------------------------------
% hpf High-pass filter
% ---------------------------------------------------------------------

hpf_wavelet_iter      = cfg_entry;
hpf_wavelet_iter.tag  = 'hpf_wavelet_iter';
hpf_wavelet_iter.name = 'Wavelet iterations';
hpf_wavelet_iter.val = {4};
%hpf_wavelet_iter.def    = @(val)nirs_get_defaults('hpf_wavelet_iter', val{:});
hpf_wavelet_iter.strtype = 'r';  
hpf_wavelet_iter.num     = [1 1]; 
hpf_wavelet_iter.help    = {'Specify wavelet iterations - default is 4 in NIRS_SPM.'}; 

hpf_wavelet = cfg_branch;
hpf_wavelet.tag     = 'hpf_wavelet';
hpf_wavelet.name    = 'Wavelet Filter';
hpf_wavelet.val     = {hpf_wavelet_iter}; 
hpf_wavelet.help    = {'Specify properties of wavelet filter'};

hpf_dct_cutoff      = cfg_entry;
hpf_dct_cutoff.tag  = 'hpf_dct_cutoff';
hpf_dct_cutoff.name = 'DCT cutoff in seconds';
hpf_dct_cutoff.val  = {128};
%hpf_dct_cutoff.def    = @(val)nirs_get_defaults('hpf_dct_cutoff', val{:});
hpf_dct_cutoff.strtype = 'r';  
hpf_dct_cutoff.num     = [1 1]; 
hpf_dct_cutoff.help    = {'Specify DCT cutoff in seconds.'}; 

hpf_dct = cfg_branch;
hpf_dct.tag     = 'hpf_dct';
hpf_dct.name    = 'DCT Filter';
hpf_dct.val     = {hpf_dct_cutoff}; 
hpf_dct.help    = {'Specify properties of Discrete Cosine Transform filter'};

hpf_none         = cfg_branch;
hpf_none.tag     = 'hpf_none';
hpf_none.name    = 'No high pass filter'; 
hpf_none.help    = {'No high pass filter.'};

nirs_hpf           = cfg_choice;
nirs_hpf.name      = 'High-pass filter';
nirs_hpf.tag       = 'nirs_hpf';
nirs_hpf.values    = {hpf_none
                      hpf_wavelet
                      hpf_dct}; 
nirs_hpf.val       = {hpf_wavelet}; 
%nirs_hpf.def    = @(val)nirs_get_defaults('nirs_hpf', val{:});
nirs_hpf.help      = {'Choose high-pass filter.'}; 

% ---------------------------------------------------------------------
% lpf Low-pass filter
% ---------------------------------------------------------------------
fwhm1      = cfg_entry;
fwhm1.tag  = 'fwhm1';
fwhm1.name = 'FWHM in seconds';
fwhm1.val = {1.5};
%fwhm1.def    = @(val)nirs_get_defaults('fwhm1', val{:});
fwhm1.strtype = 'r';  
fwhm1.num     = [1 1]; 
%fwhm1.def = @(val)nirs_get_defaults('configMC1.scalpPpties_l2', val{:});
fwhm1.help    = {'FWHM in seconds.'}; 

lpf_gauss         = cfg_branch;
lpf_gauss.tag     = 'lpf_gauss';
lpf_gauss.name    = 'Gaussian Filter';
lpf_gauss.val     = {fwhm1}; 
lpf_gauss.help    = {'Specify properties of Gaussian filter'};

lpf_none         = cfg_branch;
lpf_none.tag     = 'lpf_none';
lpf_none.name    = 'No low pass filter';
lpf_none.help    = {'No low pass filter.'};

lpf_hrf         = cfg_branch;
lpf_hrf.tag     = 'lpf_hrf';
lpf_hrf.name    = 'HRF Filter';
lpf_hrf.help    = {'HRF filter'};

nirs_lpf           = cfg_choice;
nirs_lpf.name      = 'Low-pass filter';
nirs_lpf.tag       = 'nirs_lpf';
nirs_lpf.values    = {lpf_none
                      lpf_gauss
                      lpf_hrf}; 
nirs_lpf.val       = {lpf_gauss}; 
%nirs_lpf.def    = @(val)nirs_get_defaults('nirs_lpf', val{:});
nirs_lpf.help      = {'Choose low-pass filter.'}; 


% Executable Branch
HPF_LPF      = cfg_exbranch;       
HPF_LPF.name = 'Filters';             
HPF_LPF.tag  = 'HPF_LPF'; 
HPF_LPF.val  = {NIRSmat DelPreviousData NewDirCopyNIRS nirs_hpf nirs_lpf}; 
HPF_LPF.prog = @nirs_run_HPF_LPF;  
HPF_LPF.vout = @nirs_cfg_vout_HPF_LPF; 
HPF_LPF.help = {'Filters: currently only low pass, with or without ',...
    'a downsampling factor.'};

function vout = nirs_cfg_vout_HPF_LPF(job)
vout = cfg_dep;                     
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configure input files for Monte Carlo simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

latest_mcim         = cfg_entry;
latest_mcim.tag     = 'latest_mcim';
latest_mcim.name    = 'Latest ROI';
latest_mcim.val     = {'latestROI'}; 
latest_mcim.help    = {'Latest ROI generated selected.'};

mcim_in         = cfg_files; %Select MC segmented volume for this subject 
mcim_in.name    = 'MC segmented volume'; % The displayed name
mcim_in.tag     = 'mcim_in';       %file names
mcim_in.filter = 'image';
mcim_in.ufilter = '.nii';    
mcim_in.num     = [1 1];     % Number of inputs required 
mcim_in.def    = @(val)nirs_get_defaults('configMC1.image_in', val{:});
mcim_in.help    = {'Select MC segmented volume for this subject.'}; % help text displayed

% select_mcim         = cfg_branch;
% select_mcim.tag     = 'select_mcim';
% select_mcim.name    = 'Selected image';
% select_mcim.val     = {mcim_in};
% select_mcim.help    = {'Choose image.'};

mcim_cfg           = cfg_choice;
mcim_cfg.name      = 'Image';
mcim_cfg.tag       = 'mcim_cfg';
mcim_cfg.values    = {latest_mcim mcim_in};
mcim_cfg.val       = {latest_mcim}; 
mcim_cfg.help      = {'bla'}; 

MC_CUDAchoice    = cfg_menu;
MC_CUDAchoice.name   = 'Configuration file type';
MC_CUDAchoice.tag    = 'MC_CUDAchoice';
MC_CUDAchoice.labels = {'MCX: Monte Carlo Extreme','tMCimg','Both'};
MC_CUDAchoice.values = {1,2,3};
MC_CUDAchoice.def    = @(val)nirs_get_defaults('configMC1.MC_CUDAchoice', val{:});
MC_CUDAchoice.help   = {'Choose type of configuration files to generate.'};

% est ce qu'il y a pas un pb du au fait qu'il attend un directory ???
MC_configdir         = cfg_entry;
MC_configdir.tag     = 'MC_configdir';
MC_configdir.name    = 'Monte Carlo configuration files directory';
MC_configdir.strtype = 's';
MC_configdir.num     = [1 Inf];
MC_configdir.def     = @(val)nirs_get_defaults('configMC1.MC_configdir', val{:}); 
MC_configdir.help    = {'Directory to put Monte Carlo configuration files.'};

MC_nam         = cfg_entry;
MC_nam.tag     = 'MC_nam';
MC_nam.name    = 'Monte Carlo simulation name';
MC_nam.strtype = 's';
MC_nam.num     = [1 Inf];
MC_nam.help    = {'Name of Monte Carlo simulation.'};

%--------------------------------------------------------------------------
nphotons         = cfg_entry; 
nphotons.name    = 'Number of photons'; % The displayed name
nphotons.tag     = 'nphotons';       %file names
nphotons.strtype = 'r';  
nphotons.num     = [1 1];     % Number of inputs required 
nphotons.def = @(val)nirs_get_defaults('configMC1.nphotons', val{:});
nphotons.help    = {'Input number of photons (not currently used).'}; 
             
seed         = cfg_entry; 
seed.name    = 'Random seed'; % The displayed name
seed.tag     = 'seed';       %file names
seed.strtype = 'r';  
seed.num     = [1 1];     % Number of inputs required 
seed.def = @(val)nirs_get_defaults('configMC1.seed', val{:});
seed.help    = {'Input random seed.'}; 

modulationFreq         = cfg_entry; 
modulationFreq.name    = 'Modulation Frequency'; % The displayed name
modulationFreq.tag     = 'modulationFreq';       %file names
modulationFreq.strtype = 'r';  
modulationFreq.num     = [1 1];     % Number of inputs required 
modulationFreq.def = @(val)nirs_get_defaults('configMC1.modulationFreq', val{:});
modulationFreq.help    = {'Modulation Frequency; leave at 0 for CW (continuous wave operation).'}; 

deltaT         = cfg_entry; 
deltaT.name    = 'deltaT'; % The displayed name
deltaT.tag     = 'deltaT';       %file names
deltaT.strtype = 'r';  
deltaT.num     = [1 1];     % Number of inputs required 
deltaT.def = @(val)nirs_get_defaults('configMC1.deltaT', val{:});
deltaT.help    = {'deltaT: interval between time gates in seconds.'}; 

numTimeGates         = cfg_entry; 
numTimeGates.name    = 'numTimeGates'; % The displayed name
numTimeGates.tag     = 'numTimeGates';       %file names
numTimeGates.strtype = 'r';  
numTimeGates.num     = [1 1];     % Number of inputs required 
numTimeGates.def = @(val)nirs_get_defaults('configMC1.numTimeGates', val{:});
numTimeGates.help    = {'Number of time gates; total duration is number of time gates times deltaT.'}; 

radiis         = cfg_entry; 
radiis.name    = 'Source Radii'; % The displayed name
radiis.tag     = 'radiis';       %file names
radiis.strtype = 'r';  
radiis.num     = [1 1];     % Number of inputs required 
radiis.def = @(val)nirs_get_defaults('configMC1.radiis', val{:});
radiis.help    = {'Input radius of the sources.'}; 

radiid         = cfg_entry; 
radiid.name    = 'Detector radii'; % The displayed name
radiid.tag     = 'radiid';       %file names
radiid.strtype = 'r';  
radiid.num     = [1 1];     % Number of inputs required 
radiid.def = @(val)nirs_get_defaults('configMC1.radiid', val{:});
radiid.help    = {'Input radius of the detectors.'}; 

voxelSize         = cfg_entry; 
voxelSize.name    = 'Voxel Size'; % The displayed name
voxelSize.tag     = 'voxelSize';       %file names
voxelSize.strtype = 'r';  
voxelSize.num     = [1 1];     % Number of inputs required 
voxelSize.def = @(val)nirs_get_defaults('configMC1.voxelSize', val{:});
voxelSize.help    = {'Input voxel Size.'}; 

gmPpties_l1         = cfg_entry; 
gmPpties_l1.name    = 'Gray Matter first wavelength'; % The displayed name
gmPpties_l1.tag     = 'gmPpties_l1';       %file names
gmPpties_l1.strtype = 'r';  
gmPpties_l1.num     = [1 4];     % Number of inputs required 
gmPpties_l1.def = @(val)nirs_get_defaults('configMC1.gmPpties_l1', val{:});
gmPpties_l1.help    = {'Gray matter properties (\mu_a,\mu_s, g, n) for first wavelength (default = 830 nm).'}; 

wmPpties_l1         = cfg_entry; 
wmPpties_l1.name    = 'White Matter first wavelength'; % The displayed name
wmPpties_l1.tag     = 'wmPpties_l1';       %file names
wmPpties_l1.strtype = 'r';  
wmPpties_l1.num     = [1 4];     % Number of inputs required 
wmPpties_l1.def = @(val)nirs_get_defaults('configMC1.wmPpties_l1', val{:});
wmPpties_l1.help    = {'White matter properties (\mu_a,\mu_s, g, n) for first wavelength (default = 830 nm).'}; 

csfPpties_l1         = cfg_entry; 
csfPpties_l1.name    = 'CSF first wavelength'; % The displayed name
csfPpties_l1.tag     = 'csfPpties_l1';       %file names
csfPpties_l1.strtype = 'r';  
csfPpties_l1.num     = [1 4];     % Number of inputs required 
csfPpties_l1.def = @(val)nirs_get_defaults('configMC1.csfPpties_l1', val{:});
csfPpties_l1.help    = {'CSF properties (\mu_a,\mu_s, g, n) for first wavelength (default = 830 nm).'}; 

skullPpties_l1         = cfg_entry; 
skullPpties_l1.name    = 'Skull first wavelength'; % The displayed name
skullPpties_l1.tag     = 'skullPpties_l1';       %file names
skullPpties_l1.strtype = 'r';  
skullPpties_l1.num     = [1 4];     % Number of inputs required 
skullPpties_l1.def = @(val)nirs_get_defaults('configMC1.skullPpties_l1', val{:});
skullPpties_l1.help    = {'Skull properties (\mu_a,\mu_s, g, n) for first wavelength (default = 830 nm).'}; 

scalpPpties_l1         = cfg_entry; 
scalpPpties_l1.name    = 'Scalp first wavelength'; % The displayed name
scalpPpties_l1.tag     = 'scalpPpties_l1';       %file names
scalpPpties_l1.strtype = 'r';  
scalpPpties_l1.num     = [1 4];     % Number of inputs required 
scalpPpties_l1.def = @(val)nirs_get_defaults('configMC1.scalpPpties_l1', val{:});
scalpPpties_l1.help    = {'Scalp properties (\mu_a,\mu_s, g, n) for first wavelength (default = 830 nm).'}; 

perturbationPpties_l1 = cfg_entry; 
perturbationPpties_l1.name    = 'Perturbation first wavelength'; % The displayed name
perturbationPpties_l1.tag     = 'perturbationPpties_l1';       %file names
perturbationPpties_l1.strtype = 'r';  
perturbationPpties_l1.num     = [1 4];     % Number of inputs required 
perturbationPpties_l1.def = @(val)nirs_get_defaults('configMC1.perturbationPpties_l1', val{:});
perturbationPpties_l1.help    = {['Action on grey matter only: ',...
    'Perturbation properties Delta(\mu_a,\mu_s, g, n) for first wavelength (default = 830 nm).']}; 

gmPpties_l2         = cfg_entry; 
gmPpties_l2.name    = 'Gray Matter second wavelength'; % The displayed name
gmPpties_l2.tag     = 'gmPpties_l2';       %file names
gmPpties_l2.strtype = 'r';  
gmPpties_l2.num     = [1 4];     % Number of inputs required 
gmPpties_l2.def = @(val)nirs_get_defaults('configMC1.gmPpties_l2', val{:});
gmPpties_l2.help    = {'Gray matter properties (\mu_a,\mu_s, g, n) for second wavelength (default = 690 nm).'}; 

wmPpties_l2         = cfg_entry; 
wmPpties_l2.name    = 'White Matter second wavelength'; % The displayed name
wmPpties_l2.tag     = 'wmPpties_l2';       %file names
wmPpties_l2.strtype = 'r';  
wmPpties_l2.num     = [1 4];     % Number of inputs required 
wmPpties_l2.def = @(val)nirs_get_defaults('configMC1.wmPpties_l2', val{:});
wmPpties_l2.help    = {'White matter properties (\mu_a,\mu_s, g, n) for second wavelength (default = 690 nm).'}; 

csfPpties_l2         = cfg_entry; 
csfPpties_l2.name    = 'CSF second wavelength'; % The displayed name
csfPpties_l2.tag     = 'csfPpties_l2';       %file names
csfPpties_l2.strtype = 'r';  
csfPpties_l2.num     = [1 4];     % Number of inputs required 
csfPpties_l2.def = @(val)nirs_get_defaults('configMC1.csfPpties_l2', val{:});
csfPpties_l2.help    = {'CSF properties (\mu_a,\mu_s, g, n) for second wavelength (default = 690 nm).'}; 

skullPpties_l2         = cfg_entry; 
skullPpties_l2.name    = 'Skull second wavelength'; % The displayed name
skullPpties_l2.tag     = 'skullPpties_l2';       %file names
skullPpties_l2.strtype = 'r';  
skullPpties_l2.num     = [1 4];     % Number of inputs required 
skullPpties_l2.def = @(val)nirs_get_defaults('configMC1.skullPpties_l2', val{:});
skullPpties_l2.help    = {'Skull properties (\mu_a,\mu_s, g, n) for second wavelength (default = 690 nm).'}; 

scalpPpties_l2         = cfg_entry; 
scalpPpties_l2.name    = 'Scalp second wavelength'; % The displayed name
scalpPpties_l2.tag     = 'scalpPpties_l2';       %file names
scalpPpties_l2.strtype = 'r';  
scalpPpties_l2.num     = [1 4];     % Number of inputs required 
scalpPpties_l2.def = @(val)nirs_get_defaults('configMC1.scalpPpties_l2', val{:});
scalpPpties_l2.help    = {'Scalp properties (\mu_a,\mu_s, g, n) for second wavelength (default = 690 nm).'}; 

perturbationPpties_l2 = cfg_entry; 
perturbationPpties_l2.name    = 'Perturbation second wavelength'; % The displayed name
perturbationPpties_l2.tag     = 'perturbationPpties_l2';       %file names
perturbationPpties_l2.strtype = 'r';  
perturbationPpties_l2.num     = [1 4];     % Number of inputs required 
perturbationPpties_l2.def = @(val)nirs_get_defaults('configMC1.perturbationPpties_l2', val{:});
perturbationPpties_l2.help    = {['Action on grey matter only: ',...
    'Perturbation properties Delta(\mu_a,\mu_s, g, n) for first wavelength (default = 830 nm).']}; 

MC_parameters      = cfg_branch;
MC_parameters.tag  = 'MC_parameters';
MC_parameters.name = 'Parameters';
MC_parameters.val  = {nphotons seed modulationFreq deltaT numTimeGates radiis radiid voxelSize ...
    gmPpties_l1 wmPpties_l1 csfPpties_l1 skullPpties_l1 scalpPpties_l1 perturbationPpties_l1 ...
    gmPpties_l2 wmPpties_l2 csfPpties_l2 skullPpties_l2 scalpPpties_l2 perturbationPpties_l2};
MC_parameters.help = {'Parameters'};

% Executable Branch
configMC1      = cfg_exbranch;       
configMC1.name = 'Configure Monte Carlo inputs';            
configMC1.tag  = 'configMC1'; 
configMC1.val  = {MC_nam NIRSmat NewDirCopyNIRS mcim_cfg MC_CUDAchoice MC_configdir MC_parameters};    
configMC1.prog = @nirs_run_configMC;  
configMC1.vout = @nirs_cfg_vout_configMC; 
configMC1.help = {'Generate configuration input files for Monte Carlo simulation.'};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_configMC(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration: run MC 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MCXconfigFiles         = cfg_files;  
MCXconfigFiles.name    = 'Select input files'; 
MCXconfigFiles.tag     = 'MCXconfigFiles';       
MCXconfigFiles.ufilter = '.inp';    
MCXconfigFiles.num     = [1 Inf];      
MCXconfigFiles.help    = {'Select input files (.inp for MCX).'}; 

MCX1         = cfg_branch;
MCX1.tag     = 'MCX1';
MCX1.name    = 'Monte Carlo Extreme';
MCX1.val     = {MCXconfigFiles }; % MCXconfig};
MCX1.help    = {'Run Monte Carlo Extreme simulation'};

tMCimg_configFiles         = cfg_files; %Select 
tMCimg_configFiles.name    = 'Select input files'; 
tMCimg_configFiles.tag     = 'tMCimg_configFiles';       
tMCimg_configFiles.ufilter = '.cfg';    
tMCimg_configFiles.num     = [1 Inf];     
tMCimg_configFiles.help    = {'Select input files (.cfg for tMCimg).'}; 

tMCimg1         = cfg_branch;
tMCimg1.tag     = 'tMCimg1';
tMCimg1.name    = 'tMCimg Monte Carlo Simulation';
tMCimg1.val     = {tMCimg_configFiles};
tMCimg1.help    = {'Run tMCimg Monte Carlo simulation'};
% 
MC_runCUDAchoice        = cfg_choice;
MC_runCUDAchoice.name   = 'Monte Carlo simulation method';
MC_runCUDAchoice.tag    = 'MC_runCUDAchoice';
MC_runCUDAchoice.values = {MCX1,tMCimg1};
MC_runCUDAchoice.val    = {MCX1};
MC_runCUDAchoice.help   = {['Choose method of Monte Carlo simulation. ',...
    'MCX is much faster but requires a CUDA compatible graphics card.']};

% Executable Branch
runMC1      = cfg_exbranch;      
runMC1.name = 'Run Monte Carlo simulation';           
runMC1.tag  = 'runMC1';
runMC1.val  = {MC_runCUDAchoice}; 
runMC1.prog = @nirs_run_runMC;  
runMC1.vout = @nirs_cfg_vout_runMC; 
runMC1.help = {'Run Monte Carlo simulation.'};

%make .mc2 or (.his, .2pt) file names available as a dependency
function vout = nirs_cfg_vout_runMC(job)
vout = cfg_dep;                    
vout.sname      = 'MCX or tMCimg output files';       
vout.src_output = substruct('.','outMCfiles'); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

image_in         = cfg_files;
image_in.name    = 'Anatomical segmented image';
image_in.tag     = 'image_in';
image_in.filter  = 'image'; 
image_in.num     = [1 1];
image_in.help    = {'Select the SAME image as the one used to run the MC simulations'};

out_dir         = cfg_files;
out_dir.tag     = 'out_dir';
out_dir.name    = 'Output Directory';
out_dir.help    = {'Select a directory where the output image will be written.'};
out_dir.filter = 'dir';
out_dir.ufilter = '.*';
out_dir.num     = [1 1];

out_dim      = cfg_entry;
out_dim.tag  = 'out_dim';
out_dim.name = 'Output dimension';
out_dim.val = {1};
out_dim.strtype = 'r';  
out_dim.num     = [1 3]; 
out_dim.def  = @(val)nirs_get_defaults('coregNIRS.resize1.out_dim', val{:});
out_dim.help = {['Enter output image size or let [1 1 1] if you just ',...
            'want to get isotropic voxels image.']};

out_dt      = cfg_entry;
out_dt.tag  = 'out_dt';
out_dt.name = 'Data type';
out_dt.val = {1};
out_dt.strtype = 's';
out_dt.num     = [1 Inf];
out_dt.def  = @(val)nirs_get_defaults('coregNIRS.resize1.out_dt', val{:});
out_dt.help = {['Enter output image size or let ''same'' if you just ',...
        'want to get isotropic voxels image.']};

out_autonaming      = cfg_menu;
out_autonaming.tag  = 'out_autonaming';
out_autonaming.name = 'Automatic output naming';
out_autonaming.labels = {'Yes','No'};
out_autonaming.values = {0,1};
out_autonaming.def  = @(val)nirs_get_defaults('coregNIRS.resize1.out_autonaming', val{:});
out_autonaming.help = {['Choose wheather you want to choose the name ',...
    'of the output or not. If answer is ''Yes'', please change enter name.']};

out_prefix      = cfg_entry;
out_prefix.tag  = 'out_prefix';
out_prefix.name = 'Prefix of the output image';
out_prefix.strtype = 's';
out_prefix.num     = [1 Inf];
out_prefix.help = {['You can choose to give a particular prefix to the ',...
        'output image. This prefix will be added at the left of the name ',...
        'of the image. A default name will be given otherwise.']};
    
% Executable Branch
resize1      = cfg_exbranch;
resize1.name = 'Resize image';
resize1.tag  = 'resize1';
resize1.val  = {image_in out_dir out_dim out_dt out_autonaming out_prefix};
resize1.prog = @nirs_resize;
resize1.vout = @nirs_cfg_vout_resize;
resize1.help = {'Resize the input image with respect to output size.'};


function vout = nirs_cfg_vout_resize(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration: generate sensitivity matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outMCfiles      = cfg_files;
outMCfiles.name    = 'Select MC output files';
outMCfiles.tag     = 'outMCfiles';
outMCfiles.ufilter = {'.2pt','.mc2'};    
outMCfiles.num     = [1 Inf];     
outMCfiles.help    = {'Select .mc2 or .2pt files for this subject.'}; 

% Executable Branch
makesens1      = cfg_exbranch;      
makesens1.name = 'Sensitivity Matrix';            
makesens1.tag  = 'makesens1'; 
makesens1.val  = {outMCfiles NIRSmat};
makesens1.prog = @nirs_run_generate_sensitivity_matrix;  
makesens1.vout = @nirs_cfg_vout_generate_sensitivity_matrix; 
makesens1.help = {'Generate sensitivity matrix.'};

function vout = nirs_cfg_vout_generate_sensitivity_matrix(job)
vout = cfg_dep;                     
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration: 3D reconstruction -- Tikhonov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Executable Branch
tikhonov1      = cfg_exbranch;       
tikhonov1.name = 'Tikhonov inversion';             
tikhonov1.tag  = 'tikhonov1';
tikhonov1.val  = {NIRSmat}; 
tikhonov1.prog = @nirs_run_inverse_tikhonov;  
tikhonov1.vout = @nirs_cfg_vout_inverse_tikhonov; 
tikhonov1.help = {'Invert using Tikhonov.'};

function vout = nirs_cfg_vout_inverse_tikhonov(job)
vout = cfg_dep;                     
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration: 3D reconstruction -- ReML reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sensmat         = cfg_files; %Select sensitivity matrix for this subject
% sensmat.name    = 'Sensitivity Matrix'; % The displayed name
% sensmat.tag     = 'sensmat';       %file names
% sensmat.filter  = 'nii';
% sensmat.num     = [1 1];     % Number of inputs required
% sensmat.help    = {'Select sensitivity matrix for this subject.'}; % help t

anat_segT1         = cfg_files; %Select MC segmented volume for this subject
anat_segT1.name    = 'Anatomical segmented image'; % The displayed name
anat_segT1.tag     = 'anat_segT1';       %file names
anat_segT1.filter  = 'image';
anat_segT1.ufilter = '.nii';
anat_segT1.num     = [1 1];     % Number of inputs required
anat_segT1.help    = {'Anatomical segmented image with NewSegment and MCsegment.'}; % help text displayed

% Priors
hb_relative_evolution        = cfg_menu;
hb_relative_evolution.name   = 'Relative evolution HbO/HbR';
hb_relative_evolution.tag    = 'hb_relative_evolution';
hb_relative_evolution.labels = {'Yes','No'};
hb_relative_evolution.values = {1,0};
%hb_relative_evolution.def    = @(val)nirs_get_defaults('hb_relative_evolution', val{:});
hb_relative_evolution.help   = {'Choose type of configuration files to generate.'};

priors      = cfg_branch;
priors.name = 'Priors';
priors.tag  = 'priors';
priors.val  = {anat_segT1 hb_relative_evolution};
priors.help = {'Choose priors you want to use for the reconstruction.'};

dir_in         = cfg_files;
dir_in.tag     = 'dir_in';
dir_in.name    = 'MonteCarlo output directory';
dir_in.filter = 'dir';
dir_in.ufilter = '.*';
dir_in.num     = [1 1];
dir_in.help    = {'Select the MonteCarlo simulation output directory.'};

ReML_method          = cfg_menu;
ReML_method.name      = 'ReML method';
ReML_method.tag       = 'ReML_method';
ReML_method.labels    = {'Huppert' 'SPM' 'Tikhonov'};
ReML_method.values    = {0,1,2};
ReML_method.val       = {0};
ReML_method.help      = {'Choose ReML reconstruction method.'};

% Executable Branch
ReMLreconstruct1      = cfg_exbranch;       
ReMLreconstruct1.name = '3D NIRS data ReML reconstruction';             
ReMLreconstruct1.tag  = 'ReMLreconstruct1';
ReMLreconstruct1.val  = {NIRSmat dir_in ReML_method};   
ReMLreconstruct1.prog = @nirs_run_ReMLreconstruct;  
ReMLreconstruct1.vout = @nirs_cfg_vout_ReMLreconstruct; 
ReMLreconstruct1.help = {'Run 3D NIRS data reconstruction.'};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_ReMLreconstruct(job)
vout = cfg_dep;                     
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General Linear Model Specification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------------------------------------------------------------------
% dir Directory
% ---------------------------------------------------------------------
dir         = cfg_files;
dir.tag     = 'dir';
dir.name    = 'Stats Directory';
dir.help    = {'Select a directory where the NIRS_SPM.mat files containing the specified design matrix will be written.'};
dir.filter = 'dir';
dir.ufilter = '.*';
dir.num     = [1 1];

% ---------------------------------------------------------------------
% units Units for design
% ---------------------------------------------------------------------
units         = cfg_menu;
units.tag     = 'units';
units.name    = 'Units for design';
units.help    = {'The onsets of events or blocks can be specified in either scans or seconds.'};
units.labels = {
                'Scans'
                'Seconds'
                };
units.values  = {
                    0
                    1
                    };
units.val = {1};
%units.def = @(val)nirs_get_defaults('model_specify.units', val{:});

time_res      = cfg_entry;
time_res.tag  = 'time_res';
time_res.name = 'Time resolution';
time_res.val = {1};
time_res.strtype = 'r';  
time_res.num     = [1 1]; 
%time_res.def = @(val)nirs_get_defaults('model_specify.time_res', val{:});
time_res.help    = {'Time resolution for onsets will be given by NIRS sampling rate divided by this factor  - value is 10 in NIRS_SPM.'}; 

input_onsets         = cfg_files;  
input_onsets.name    = 'Select onset files'; % The displayed name
input_onsets.tag     = 'input_onsets';       
input_onsets.filter  = 'mat';    
input_onsets.val{1}  = {''};
input_onsets.num     = [0 Inf];     % Number of inputs required 
input_onsets.help    = {'Select onset files for each session of this subject.'}; % help text displayed

% ---------------------------------------------------------------------
% multi_reg Multiple regressors
% ---------------------------------------------------------------------
multi_reg         = cfg_files;
multi_reg.tag     = 'multi_reg';
multi_reg.name    = 'Multiple regressors';
multi_reg.val{1} = {''};
multi_reg.help    = {
                     'Select the *.mat/*.txt file containing details of your multiple regressors. '
                     ''
                     'If you have multiple regressors eg. realignment parameters, then entering the details a regressor at a time is very inefficient. This option can be used to load all the required information in one go. '
                     ''
                     'You will first need to create a *.mat file containing a matrix R or a *.txt file containing the regressors. Each column of R will contain a different regressor. When SPM creates the design matrix the regressors will be named R1, R2, R3, ..etc.'
}';
multi_reg.filter = 'mat';
multi_reg.ufilter = '.*';
multi_reg.num     = [0 Inf];


subj         = cfg_branch;
subj.tag     = 'subj';
subj.name    = 'Subject';
subj.val     = {input_onsets multi_reg};
subj.help    = {};



% ---------------------------------------------------------------------
% Noise method
% ---------------------------------------------------------------------
nirs_noise         = cfg_menu;
nirs_noise.tag     = 'nirs_noise';
nirs_noise.name    = 'Noise method';
nirs_noise.help    = {'Choose method for noise treatment.'}';
nirs_noise.labels  = {
                    'precoloring'
                    'prewhitening'
                    };
nirs_noise.values  = {
                    0
                    1
                    };
%nirs_noise.val = {0};
nirs_noise.def = @(val)nirs_get_defaults(...
    'model_specify.wls_bglm_specify.wls_or_bglm.NIRS_SPM.nirs_noise', val{:});    

% ---------------------------------------------------------------------
% derivs Model derivatives
% ---------------------------------------------------------------------
derivs         = cfg_menu;
derivs.tag     = 'derivs';
derivs.name    = 'Model derivatives';
derivs.help    = {'Model HRF Derivatives. The canonical HRF combined with time and dispersion derivatives comprise an ''informed'' basis set, as the shape of the canonical response conforms to the hemodynamic response that is commonly observed. The incorporation of the derivate terms allow for variations in subject-to-subject and voxel-to-voxel responses. The time derivative allows the peak response to vary by plus or minus a second and the dispersion derivative allows the width of the response to vary. The informed basis set requires an SPM{F} for inference. T-contrasts over just the canonical are perfectly valid but assume constant delay/dispersion. The informed basis set compares favourably with eg. FIR bases on many data sets. '};
derivs.labels = {
                 'No derivatives'
                 'Time derivatives'
                 'Time and Dispersion derivatives'
}';
derivs.values = {[0 0] [1 0] [1 1]};
derivs.val = {[0 0]};
%derivs.def = @(val)nirs_get_defaults('model_specify.derivs', val{:});

% ---------------------------------------------------------------------
% volt Model Interactions (Volterra)
% ---------------------------------------------------------------------
volt         = cfg_menu;
volt.tag     = 'volt';
volt.name    = 'Model Interactions (Volterra)';
volt.help    = {
                'Generalized convolution of inputs (U) with basis set (bf).'
                ''
                'For first order expansions the causes are simply convolved (e.g. stick functions) in U.u by the basis functions in bf to create a design matrix X.  For second order expansions new entries appear in ind, bf and name that correspond to the interaction among the orginal causes. The basis functions for these efects are two dimensional and are used to assemble the second order kernel. Second order effects are computed for only the first column of U.u.'
                'Interactions or response modulations can enter at two levels.  Firstly the stick function itself can be modulated by some parametric variate (this can be time or some trial-specific variate like reaction time) modeling the interaction between the trial and the variate or, secondly interactions among the trials themselves can be modeled using a Volterra series formulation that accommodates interactions over time (and therefore within and between trial types).'
}';
volt.labels = {
               'Do not model Interactions'
               'Model Interactions (2nd Volterra)'
               'Model 3rd Volterra'
}';
volt.values = {1 2 3};
volt.val = {2};
%volt.def = @(val)nirs_get_defaults('model_specify.volt', val{:}); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LIOM General Linear Model Specification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir1         = cfg_entry; 
dir1.name    = 'Stats Directory';
dir1.tag     = 'dir1';       
dir1.strtype = 's';
dir1.num     = [1 Inf];     
dir1.val{1}  = 'Stat';
%dir1.def    = @(val)nirs_get_defaults('readNIRS.boxy1.cf1.sizebloc', val{:}); 
dir1.help    = {'Enter a subdirectory name where the NIRS_SPM.mat files '
        'containing the specified design matrix will be written.'}'; 

LiomDeleteLarge      = cfg_menu;
LiomDeleteLarge.tag  = 'LiomDeleteLarge';
LiomDeleteLarge.name = 'Delete large files';
LiomDeleteLarge.labels = {'Yes','No'};
LiomDeleteLarge.values = {1,0};
LiomDeleteLarge.def = @(val)nirs_get_defaults(...
    'model_specify.wls_bglm_specify.LiomDeleteLarge', val{:}); 
LiomDeleteLarge.help = {'Delete large files (.nir and NIRS.mat) after each estimation.'};

GenerateHbT      = cfg_menu;
GenerateHbT.tag  = 'GenerateHbT';
GenerateHbT.name = 'Generate HbT';
GenerateHbT.labels = {'Yes','No'};
GenerateHbT.values = {1,0};
GenerateHbT.def = @(val)nirs_get_defaults(...
    'model_specify.wls_bglm_specify.GenerateHbT', val{:}); 
GenerateHbT.help = {'Generate HbT.'};

flag_window      = cfg_menu;
flag_window.tag  = 'flag_window';
flag_window.name = 'Show Design Matrix';
flag_window.labels = {'Yes','No'};
flag_window.values = {1,0};
flag_window.def = @(val)nirs_get_defaults(...
    'model_specify.wls_bglm_specify.flag_window', val{:}); 
flag_window.help = {'Show design matrix.'};

WLS_J0         = cfg_entry; 
WLS_J0.name    = 'Wavelet depth J0';
WLS_J0.tag     = 'WLS_J0';       
WLS_J0.strtype = 'r';
WLS_J0.num     = [1 1];     
WLS_J0.def     = @(val)nirs_get_defaults(...
    'model_specify.wls_bglm_specify.wls_or_bglm.WLS.WLS_J0', val{:}); 
WLS_J0.help    = {'Enter wavelet depth J0.'};

WLS_L0         = cfg_entry; 
WLS_L0.name    = 'Wavelet depth L0';
WLS_L0.tag     = 'WLS_L0';       
WLS_L0.strtype = 'r';
WLS_L0.num     = [1 1];     
WLS_L0.def     = @(val)nirs_get_defaults(...
    'model_specify.wls_bglm_specify.wls_or_bglm.WLS.WLS_L0', val{:}); 
WLS_L0.help    = {'Enter wavelet depth L0.'};

WLS_threshold_drift         = cfg_entry; 
WLS_threshold_drift.name    = 'Wavelet correlation threshold for drifts';
WLS_threshold_drift.tag     = 'WLS_threshold_drift';       
WLS_threshold_drift.strtype = 'r';
WLS_threshold_drift.num     = [1 1];     
WLS_threshold_drift.def     = @(val)nirs_get_defaults(...
    'model_specify.wls_bglm_specify.wls_or_bglm.WLS.WLS_threshold_drift', val{:}); 
WLS_threshold_drift.help    = {'Enter wavelet correlation threshold for drifts.'};

WLS         = cfg_branch;
WLS.tag     = 'WLS';
WLS.name    = 'Wavelet least-squares';
WLS.val     = {WLS_J0 WLS_threshold_drift WLS_L0}; 
WLS.help    = {'Specify options for wavelet least-squares method.'};

BGLM_fmax         = cfg_entry; 
BGLM_fmax.name    = 'Maximum frequency for drifts';
BGLM_fmax.tag     = 'BGLM_fmax';       
BGLM_fmax.strtype = 'r';
BGLM_fmax.num     = [1 1];     
BGLM_fmax.def     = @(val)nirs_get_defaults(...
    'model_specify.wls_bglm_specify.wls_or_bglm.BGLM.BGLM_fmax', val{:}); 
BGLM_fmax.help    = {'Enter maximum frequency for drifts in Hz.'};

BGLM_degre         = cfg_entry; 
BGLM_degre.name    = 'Polynomial degree for drifts';
BGLM_degre.tag     = 'BGLM_degre';       
BGLM_degre.strtype = 'r';
BGLM_degre.num     = [1 1];     
BGLM_degre.def     = @(val)nirs_get_defaults(...
    'model_specify.wls_bglm_specify.wls_or_bglm.BGLM.BGLM_degre', val{:}); 
BGLM_degre.help    = {'Enter polynomial degree for drifts.'};

BGLM_threshold_drift         = cfg_entry; 
BGLM_threshold_drift.name    = 'Threshold for drifts';
BGLM_threshold_drift.tag     = 'BGLM_threshold_drift';       
BGLM_threshold_drift.strtype = 'r';
BGLM_threshold_drift.num     = [1 1];     
BGLM_threshold_drift.def     = @(val)nirs_get_defaults(...
    'model_specify.wls_bglm_specify.wls_or_bglm.BGLM.BGLM_threshold_drift', val{:}); 
BGLM_threshold_drift.help    = {'Enter correlation threshold for drifts.'};

BGLM         = cfg_branch;
BGLM.tag     = 'BGLM';
BGLM.name    = 'Bayesian GLM';
BGLM.val     = {BGLM_fmax BGLM_degre BGLM_threshold_drift}; 
BGLM.help    = {'Specify options for Bayesian GLM method.'};

NIRS_SPM         = cfg_branch;
NIRS_SPM.tag     = 'NIRS_SPM';
NIRS_SPM.name    = 'NIRS_SPM MDL';
NIRS_SPM.val     = {nirs_noise nirs_hpf nirs_lpf}; 
NIRS_SPM.help    = {'Specify options for NIRS_SPM minimum description length(MDL).'};

wls_or_bglm      = cfg_choice;
wls_or_bglm.tag  = 'wls_or_bglm';
wls_or_bglm.name = 'WLS, BGLM, NIRS_SPM';
%wls_or_bglm.labels = {'WLS','BGLM', 'NIRS_SPM'};
%wls_or_bglm.values = {1,2,3};
wls_or_bglm.values = {WLS,BGLM,NIRS_SPM};
wls_or_bglm.val  = {NIRS_SPM};
%wls_or_bglm.def = @(val)nirs_get_defaults('model_specify.wls_bglm_specify.wls_or_bglm', val{:}); 
wls_or_bglm.help = {'Choose which GLM method to use:'
            'WLS: wavelet least square'
            'BGLM: Bayesian general linear model'
            'NIRS_SPM: Ye et al methods (MDL), with either precoloring or prewhitening.'}';
    
GLM_include_cardiac    = cfg_menu;
GLM_include_cardiac.name   = 'Include cardiac regressor';
GLM_include_cardiac.tag    = 'GLM_include_cardiac';
GLM_include_cardiac.labels = {'Yes','No'};
GLM_include_cardiac.values = {1,0};
GLM_include_cardiac.def    = @(val)nirs_get_defaults(...
    'model_specify.wls_bglm_specify.GLM_include_cardiac', val{:});
GLM_include_cardiac.help   = {'Include cardiac regressor if available.'};

GLM_include_Mayer    = cfg_menu;
GLM_include_Mayer.name   = 'Include Mayer wave regressor';
GLM_include_Mayer.tag    = 'GLM_include_Mayer';
GLM_include_Mayer.labels = {'Yes','No'};
GLM_include_Mayer.values = {1,0};
GLM_include_Mayer.def    = @(val)nirs_get_defaults(...
    'model_specify.wls_bglm_specify.GLM_include_Mayer', val{:});
GLM_include_Mayer.help   = {'Include Mayer wave regressor if available.'};

channel_pca      = cfg_menu;
channel_pca.tag  = 'channel_pca';
channel_pca.name = 'Spatial Principal Component Removal';
channel_pca.labels = {'Yes','No'};
channel_pca.values = {1,0};
channel_pca.def = @(val)nirs_get_defaults('model_specify.wls_bglm_specify.channel_pca', val{:}); 
channel_pca.help = {'Choose whether to do a channel PCA removal: '
            'Principal component analysis and removing the largest eigenvalue.'}';

lpf_butter_freq         = cfg_entry; 
lpf_butter_freq.name    = 'Cutoff frequency for LPF';
lpf_butter_freq.tag     = 'lpf_butter_freq';       
lpf_butter_freq.strtype = 'r';
lpf_butter_freq.num     = [1 1];     
lpf_butter_freq.def     = @(val)nirs_get_defaults(...
    'model_specify.wls_bglm_specify.lpf_butter.lpf_butter_On.lpf_butter_freq', val{:}); 
lpf_butter_freq.help    = {'Enter cutoff frequency in Hz for Butterworth LPF.'};

lpf_butter_On         = cfg_branch;
lpf_butter_On.tag     = 'lpf_butter_On';
lpf_butter_On.name    = 'Butterworth LP filter';
lpf_butter_On.val     = {lpf_butter_freq}; 
lpf_butter_On.help    = {'Butterworth low-pass filter.'};

lpf_butter_Off         = cfg_branch;
lpf_butter_Off.tag     = 'lpf_butter_Off';
lpf_butter_Off.name    = 'LP filter off';
lpf_butter_Off.val     = {}; 
lpf_butter_Off.help    = {'Low pass filter turned off.'};

lpf_butter      = cfg_choice;
lpf_butter.tag  = 'lpf_butter';
lpf_butter.name = 'Butterworth Low Pass Filter';
%lpf_butter.labels = {'Yes','No'};
lpf_butter.values = {lpf_butter_On lpf_butter_Off};
lpf_butter.val = {lpf_butter_Off};
%lpf_butter.def = @(val)nirs_get_defaults('model_specify.wls_bglm_specify.lpf_butter', val{:}); 
lpf_butter.help = {'Choose whether to include a Butterworth Low Pass Filter.'
        'Parameters are: order 3.'}';

hpf_butter_freq         = cfg_entry; 
hpf_butter_freq.name    = 'Cutoff frequency for HPF';
hpf_butter_freq.tag     = 'hpf_butter_freq';       
hpf_butter_freq.strtype = 'r';
hpf_butter_freq.num     = [1 1];     
hpf_butter_freq.def     = @(val)nirs_get_defaults(...
    'model_specify.wls_bglm_specify.hpf_butter.hpf_butter_On.hpf_butter_freq', val{:}); 
hpf_butter_freq.help    = {'Enter cutoff frequency in Hz for Butterworth HPF.'};

hpf_butter_On         = cfg_branch;
hpf_butter_On.tag     = 'hpf_butter_On';
hpf_butter_On.name    = 'Butterworth HP filter';
hpf_butter_On.val     = {hpf_butter_freq}; 
hpf_butter_On.help    = {'Butterworth high-pass filter.'};

hpf_butter_Off         = cfg_branch;
hpf_butter_Off.tag     = 'hpf_butter_Off';
hpf_butter_Off.name    = 'HP filter off';
hpf_butter_Off.val     = {}; 
hpf_butter_Off.help    = {'High pass filter turned off.'};

generate_trRV      = cfg_menu;
generate_trRV.tag  = 'generate_trRV';
generate_trRV.name = 'Generate TrRV';
generate_trRV.labels = {'Yes','No'};
generate_trRV.values = {1,0};
generate_trRV.val = {1};
generate_trRV.help = {'Generate TrRV and TrRVRV - needed for interpolated maps.'
    'Careful! Note that TrRV is required for the NIRS_SPM method, to calculate t-stats.'}';

filter_design_matrix      = cfg_menu;
filter_design_matrix.tag  = 'filter_design_matrix';
filter_design_matrix.name = 'Filter the design matrix';
filter_design_matrix.labels = {'Yes','No'};
filter_design_matrix.values = {1,0};
filter_design_matrix.val = {0};
filter_design_matrix.help = {'Currently under testing. Potential problem:'
    'introduces long range correlations in the design matrix that falsify'
    'the calculation of the nubmer of degrees of freedom, and thus the covariance.'}';

hpf_butter      = cfg_choice;
hpf_butter.tag  = 'hpf_butter';
hpf_butter.name = 'Butterworth High Pass Filter';
%hpf_butter.labels = {'Yes','No'};
hpf_butter.values = {hpf_butter_On hpf_butter_Off};
hpf_butter.val = {hpf_butter_On};
%hpf_butter.def = @(val)nirs_get_defaults('model_specify.wls_bglm_specify.hpf_butter', val{:}); 
hpf_butter.help = {'Choose whether to include a Butterworth High Pass Filter.'
        'Parameters are: order 3.'}';
    
% Executable Branch
wls_bglm_specify      = cfg_exbranch;       
wls_bglm_specify.name = 'LIOM GLM Specification';            
wls_bglm_specify.tag  = 'wls_bglm_specify'; 
wls_bglm_specify.val  = {NIRSmat dir1 subj units time_res derivs ...
    volt GLM_include_cardiac GLM_include_Mayer GenerateHbT flag_window ...
    channel_pca hpf_butter lpf_butter generate_trRV filter_design_matrix ...
     wls_or_bglm LiomDeleteLarge}; 
wls_bglm_specify.prog = @nirs_run_wls_bglm_specify;  
wls_bglm_specify.vout = @nirs_cfg_vout_wls_bglm_specify; 
wls_bglm_specify.help = {'Specify LIOM General Linear Model.'};

function vout = nirs_cfg_vout_wls_bglm_specify(job)
vout = cfg_dep;                     
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NIRS_SPM General Linear Model Specification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Executable Branch
NIRS_SPM_specify      = cfg_exbranch;       
NIRS_SPM_specify.name = 'NIRS_SPM GLM Specification';           
NIRS_SPM_specify.tag  = 'NIRS_SPM_specify'; 
NIRS_SPM_specify.val  = {NIRSmat dir subj units time_res derivs ...
    volt nirs_noise nirs_hpf nirs_lpf}; 
NIRS_SPM_specify.prog = @nirs_run_NIRS_SPM_specify;  
NIRS_SPM_specify.vout = @nirs_cfg_vout_NIRS_SPM_specify; 
NIRS_SPM_specify.help = {'NIRS_SPM GLM Specification.'};

function vout = nirs_cfg_vout_NIRS_SPM_specify(job)
vout = cfg_dep;                     
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NIRS_SPM General Linear Model Specification - NEW version - for batch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Executable Branch
NIRS_SPM_specify_batch      = cfg_exbranch;       
NIRS_SPM_specify_batch.name = 'NIRS_SPM GLM Specification (NEW, for batch)';           
NIRS_SPM_specify_batch.tag  = 'NIRS_SPM_specify_batch'; 
NIRS_SPM_specify_batch.val  = {NIRSmat dir1 subj units time_res derivs ...
    volt nirs_noise nirs_hpf nirs_lpf}; 
NIRS_SPM_specify_batch.prog = @nirs_run_NIRS_SPM_specify_batch;  
NIRS_SPM_specify_batch.vout = @nirs_cfg_vout_NIRS_SPM_specify_batch; 
NIRS_SPM_specify_batch.help = {'NIRS_SPM GLM Specification.'};

function vout = nirs_cfg_vout_NIRS_SPM_specify_batch(job)
vout = cfg_dep;                     
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LIOM General Linear Model Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NIRS_SPM_which_GLM      = cfg_menu;
NIRS_SPM_which_GLM.tag  = 'NIRS_SPM_which_GLM';
NIRS_SPM_which_GLM.name = 'Which GLM to estimate';
NIRS_SPM_which_GLM.labels = {'First','All', 'Last'};
NIRS_SPM_which_GLM.values = {1,2,3};
NIRS_SPM_which_GLM.val = {1};
NIRS_SPM_which_GLM.help = {'Choose which GLM (if more than one available) to estimate.'};

% Executable Branch
wls_bglm_estimate      = cfg_exbranch;       
wls_bglm_estimate.name = 'LIOM GLM Estimation';             
wls_bglm_estimate.tag  = 'wls_bglm_estimate';
wls_bglm_estimate.val  = {NIRSmat NIRS_SPM_which_GLM}; 
wls_bglm_estimate.prog = @nirs_run_wls_bglm_estimate;  
wls_bglm_estimate.vout = @nirs_cfg_vout_wls_bglm_estimate; 
wls_bglm_estimate.help = {'LIOM GLM Estimation: WLS (wavelet least square)'
            'and Bayesian GLM.'}';

function vout = nirs_cfg_vout_wls_bglm_estimate(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NIRS_SPM General Linear Model Estimation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dmx_files         = cfg_files;  
Dmx_files.name    = 'Select design matrix files'; 
Dmx_files.tag     = 'Dmx_files';       
Dmx_files.filter  = 'mat';    
Dmx_files.num     = [1 Inf];     
Dmx_files.help    = {'Select design matrix files to estimate.'};

% Executable Branch
NIRS_SPM_estimate      = cfg_exbranch;       
NIRS_SPM_estimate.name = 'NIRS_SPM GLM Estimation';            
NIRS_SPM_estimate.tag  = 'NIRS_SPM_estimate';
NIRS_SPM_estimate.val  = {Dmx_files}; 
NIRS_SPM_estimate.prog = @nirs_run_NIRS_SPM_estimate;  
NIRS_SPM_estimate.vout = @nirs_cfg_vout_NIRS_SPM_estimate; 
NIRS_SPM_estimate.help = {'NIRS_SPM GLM Estimation.'};

function vout = nirs_cfg_vout_NIRS_SPM_estimate(job)
vout = cfg_dep;                     
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NIRS_SPM HPF and LPF 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Executable Branch
NIRS_SPM_HPF_LPF      = cfg_exbranch;       
NIRS_SPM_HPF_LPF.name = 'Filters';             
NIRS_SPM_HPF_LPF.tag  = 'NIRS_SPM_HPF_LPF'; 
NIRS_SPM_HPF_LPF.val  = {NIRSmat nirs_hpf nirs_lpf}; 
NIRS_SPM_HPF_LPF.prog = @nirs_run_NIRS_SPM_HPF_LPF;  
NIRS_SPM_HPF_LPF.vout = @nirs_cfg_vout_NIRS_SPM_HPF_LPF; 
NIRS_SPM_HPF_LPF.help = {'Filters from NIRS_SPM. Note that design matrix'
    'is required for wavelet method.'}';

function vout = nirs_cfg_vout_NIRS_SPM_HPF_LPF(job)
vout = cfg_dep;                     
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NIRS_SPM General Linear Model Estimation - NEW version, for batch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Executable Branch
NIRS_SPM_estimate_batch      = cfg_exbranch;       
NIRS_SPM_estimate_batch.name = 'NIRS_SPM GLM Estimation (NEW, batch)';            
NIRS_SPM_estimate_batch.tag  = 'NIRS_SPM_estimate_batch';
NIRS_SPM_estimate_batch.val  = {NIRSmat NIRS_SPM_which_GLM}; 
NIRS_SPM_estimate_batch.prog = @nirs_run_NIRS_SPM_estimate_batch;  
NIRS_SPM_estimate_batch.vout = @nirs_cfg_vout_NIRS_SPM_estimate_batch; 
NIRS_SPM_estimate_batch.help = {'NIRS_SPM GLM Estimation.'};

function vout = nirs_cfg_vout_NIRS_SPM_estimate_batch(job)
vout = cfg_dep;                     
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NIRS_SPM Contrast calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%File for channel coregistration 'preproc_info' generated by NIRS_SPM
NIRS_SPM_Coregistration_Channels         = cfg_files;  
NIRS_SPM_Coregistration_Channels.name    = 'Select file of coregistration info'; 
NIRS_SPM_Coregistration_Channels.tag     = 'NIRS_SPM_Coregistration_Channels';       
NIRS_SPM_Coregistration_Channels.filter  = 'mat';    
NIRS_SPM_Coregistration_Channels.num     = [1 1];    
NIRS_SPM_Coregistration_Channels.help    = {'Select file of channel '
    'coregistration ''preproc_info'' generated by NIRS_SPM.'}'; 

%Select view
view         = cfg_entry; 
view.name    = 'View'; 
view.tag     = 'view';    
view.strtype = 'r'; 
view.num     = [1 Inf];    
view.help    = {['Enter view.  ',...
    '1: ventral  ',...
    '2: dorsal  ',...
    '3: right  ',...
    '4: left  ',...
    '5: frontal  ',...
    '6: occipital']}; % help text displayed

%Contrast name
contrast_name         = cfg_entry; 
contrast_name.name    = 'Contrast name';
contrast_name.tag     = 'contrast_name';       
contrast_name.strtype = 's';
contrast_name.num     = [1 Inf];     
contrast_name.help    = {'Contrast name'}; 

%Contrast vector
contrast_c         = cfg_entry; 
contrast_c.name    = 'Contrast vector';
contrast_c.tag     = 'contrast_c';       
contrast_c.strtype = 'r';
contrast_c.num     = [1 Inf];     
contrast_c.help    = {'Contrast vector'}; 

contrast_data         = cfg_branch;
contrast_data.tag     = 'contrast_data';
contrast_data.name    = 'Contrasts';
contrast_data.val     = {contrast_name contrast_c}; %contrast_type
contrast_data.help    = {'Specify contrasts.'};

contrast_struct         = cfg_repeat;
contrast_struct.tag     = 'contrast_struct';
contrast_struct.name    = 'Contrasts';
contrast_struct.help    = {'Specify contrasts'};
contrast_struct.values  = {contrast_data};
contrast_struct.num     = [1 Inf];

% Executable Branch
NIRS_SPM_contrast      = cfg_exbranch;      
NIRS_SPM_contrast.name = 'NIRS_SPM Contrast Calculations';            
NIRS_SPM_contrast.tag  = 'NIRS_SPM_contrast';
NIRS_SPM_contrast.val  = {Dmx_files NIRS_SPM_Coregistration_Channels ...
                view contrast_struct}; 
NIRS_SPM_contrast.prog = @nirs_run_NIRS_SPM_contrast;  
NIRS_SPM_contrast.vout = @nirs_cfg_vout_NIRS_SPM_contrast; 
NIRS_SPM_contrast.help = {'NIRS_SPM Contrast Calculations.'};

function vout = nirs_cfg_vout_NIRS_SPM_contrast(job)
vout = cfg_dep;                     
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Liom Contrast calculations - based on tube formula and code by NIRS_SPM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

liom_contrast_struct         = cfg_repeat;
liom_contrast_struct.tag     = 'liom_contrast_struct';
liom_contrast_struct.name    = 'Contrasts';
liom_contrast_struct.help    = {'Specify contrasts'};
liom_contrast_struct.values  = {contrast_data};
liom_contrast_struct.num     = [0 Inf];

contrast_p_value         = cfg_entry; 
contrast_p_value.name    = 'Contrast uncorrected p_value';
contrast_p_value.tag     = 'contrast_p_value';       
contrast_p_value.strtype = 'r';
contrast_p_value.num     = [1 1];
contrast_p_value.val     = {0.05};
contrast_p_value.help    = {'Contrast uncorrected p_value'}; 

contrast_figures      = cfg_menu;
contrast_figures.tag  = 'contrast_figures';
contrast_figures.name = 'Generate figures';
contrast_figures.labels = {'No','Both .fig and .tiff','Only .fig','Only .tiff'};
contrast_figures.values = {0,1,2,3};
contrast_figures.val = {0};
contrast_figures.help = {'Generate contrast figures. '
    'Note .fig colorbar is incorrect - it is not saved properly by Matlab.'
    'Use .tiff to view colorbar for t-stat.'}';

figures_visible      = cfg_menu;
figures_visible.tag  = 'figures_visible';
figures_visible.name = 'Make figures visible';
figures_visible.labels = {'Yes','No'};
figures_visible.values = {1,0};
figures_visible.val = {0};
figures_visible.help = {'Make figures visible during processing.'}';

colorbar_max         = cfg_entry; 
colorbar_max.name    = 'Colorbar maximum value';
colorbar_max.tag     = 'colorbar_max';       
colorbar_max.strtype = 'r';
colorbar_max.num     = [1 1];
colorbar_max.val     = {5};
colorbar_max.help    = {'Enter maximum value for colorbar'}; 

colorbar_min         = cfg_entry; 
colorbar_min.name    = 'Colorbar minimum value';
colorbar_min.tag     = 'colorbar_min';       
colorbar_min.strtype = 'r';
colorbar_min.num     = [1 1];
colorbar_min.val     = {2};
colorbar_min.help    = {'Enter minimum value for colorbar'}; 

colorbar_override      = cfg_branch;
colorbar_override.name      = 'Override colorbar';
colorbar_override.tag       = 'colorbar_override';
colorbar_override.val       = {colorbar_min colorbar_max}; 
colorbar_override.help      = {'Override colorbar.'};

colorbar_default      = cfg_branch;
colorbar_default.name      = 'Default colorbar';
colorbar_default.tag       = 'colorbar_default';
colorbar_default.val       = {}; 
colorbar_default.help      = {'Default colorbar.'};

override_colorbar           = cfg_choice;
override_colorbar.name      = 'Override colorbar';
override_colorbar.tag       = 'override_colorbar';
override_colorbar.values    = {colorbar_default colorbar_override};
override_colorbar.val       = {colorbar_default}; 
override_colorbar.help      = {'Override default treatment of colorbar.'
    'User can then specify maximum and minimum values for the colorbar.'}';

GenerateInverted      = cfg_menu;
GenerateInverted.tag  = 'GenerateInverted';
GenerateInverted.name = 'Generate Inverted Responses';
GenerateInverted.labels = {'Yes','No'};
GenerateInverted.values = {1,0};
GenerateInverted.val = {1};
GenerateInverted.help = {'Generate contrasts for inverted responses.'};

GroupFiguresIntoSubplots      = cfg_menu;
GroupFiguresIntoSubplots.tag  = 'GroupFiguresIntoSubplots';
GroupFiguresIntoSubplots.name = 'Group Figures Into Subplots';
GroupFiguresIntoSubplots.labels = {'Yes','No'};
GroupFiguresIntoSubplots.values = {1,0};
GroupFiguresIntoSubplots.val = {1};
GroupFiguresIntoSubplots.help = {'Group Figures Into Subplots.'};

% Executable Branch
liom_contrast      = cfg_exbranch;      
liom_contrast.name = 'Liom Contrast Calculations';            
liom_contrast.tag  = 'liom_contrast';
liom_contrast.val  = {NIRSmat view liom_contrast_struct GenerateInverted contrast_p_value ...
    contrast_figures override_colorbar figures_visible GroupFiguresIntoSubplots TopoData}; 
liom_contrast.prog = @nirs_run_liom_contrast;  
liom_contrast.vout = @nirs_cfg_vout_liom_contrast; 
liom_contrast.help = {'Liom Contrast Calculations.'};

function vout = nirs_cfg_vout_liom_contrast(job)
vout = cfg_dep;                     
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NIRS_SPM Group Level Model Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Contrast_files         = cfg_files;  
Contrast_files.name    = 'Select estimated contrasts files'; 
Contrast_files.tag     = 'Contrast_files';       
Contrast_files.filter  = 'mat';    
Contrast_files.num     = [1 Inf];     
Contrast_files.help    = {'Select estimated constrast files for this '
            'group. Select all desired files, and code will try to group '
            'them by view type, contrast type, and by chromophore. '
            'Please see code if any doubt.'}'; 

% Executable Branch
NIRS_SPM_group      = cfg_exbranch;       
NIRS_SPM_group.name = 'NIRS_SPM Group Model Estimation';             
NIRS_SPM_group.tag  = 'NIRS_SPM_group'; 
NIRS_SPM_group.val  = {Contrast_files}; 
NIRS_SPM_group.prog = @nirs_run_NIRS_SPM_group;  
NIRS_SPM_group.vout = @nirs_cfg_vout_NIRS_SPM_group; 
NIRS_SPM_group.help = {'NIRS_SPM Group level model estimation.'};

function vout = nirs_cfg_vout_NIRS_SPM_group(job)
vout = cfg_dep;                     
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Liom Group Level Model Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dir_ga         = cfg_entry; 
% dir_ga.name    = 'Group Analysis';
% dir_ga.tag     = 'dir_ga';       
% dir_ga.strtype = 's';
% dir_ga.num     = [1 Inf];     
% dir_ga.val{1}  = 'Stat';
% dir_ga.help    = {'Enter a subdirectory name where the NIRS_SPM.mat files '
%         'containing the specified design matrix will be written.'}'; 

session_number         = cfg_entry; 
session_number.name    = 'Session number';
session_number.tag     = 'session_number';       
session_number.strtype = 'r';
session_number.num     = [1 1];     
session_number.help    = {'Enter the number of the session you want to analyse. Only one session can be analysed at a time.'}; 

FFX_or_RFX = cfg_menu;
FFX_or_RFX.tag  = 'FFX_or_RFX';
FFX_or_RFX.name = 'Fixed or random effects';
FFX_or_RFX.labels = {'FFX','RFX'};
FFX_or_RFX.values = {1,0};
FFX_or_RFX.val = {1};
FFX_or_RFX.help = {'Use fixed effects (FFX) for group of sessions (intra-subject) '
    'Use random effects (RFX) for group of subjects (inter-subject)'
    'FFX amounts to setting the between session variance to zero.'
    'Several subjects can be specified for FFX; they will each be treated separately.'
    'RFX assumes there is only one session per subject.'
    'RFX can also be used for one subject with multiple sessions, '
    'to take into account the variance between sessions.'}';

% Executable Branch
liom_group      = cfg_exbranch;       
liom_group.name = 'Liom Group Model Estimation';             
liom_group.tag  = 'liom_group'; 
liom_group.val  = {NIRSmat FFX_or_RFX contrast_figures contrast_p_value ...
        GenerateInverted override_colorbar figures_visible GroupFiguresIntoSubplots}; 
liom_group.prog = @nirs_run_liom_group;  
liom_group.vout = @nirs_cfg_vout_liom_group; 
liom_group.help = {'Liom Group level model estimation.'};

function vout = nirs_cfg_vout_liom_group(job)
vout = cfg_dep;                     
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyze GLMs - loop over subjects and jobs, to plot simple t contrasts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ROCLoopJob         = cfg_files;  
ROCLoopJob.name    = 'Select job(s) to loop over'; 
ROCLoopJob.tag     = 'ROCLoopJob';       
ROCLoopJob.ufilter = '.mat';   
ROCLoopJob.num     = [1 Inf];    
ROCLoopJob.help    = {'Select .mat-format previously specified '
                        'and saved job(s) to loop over.'}'; 

% Executable Branch
AnalyzeGLM      = cfg_exbranch;       
AnalyzeGLM.name = 'Analyze GLMs';            
AnalyzeGLM.tag  = 'AnalyzeGLM';
AnalyzeGLM.val  = {NIRSmat ROCLoopJob}; 
AnalyzeGLM.prog = @nirs_run_AnalyzeGLM;  
AnalyzeGLM.vout = @nirs_cfg_vout_AnalyzeGLM; 
AnalyzeGLM.help = {'This module performs a large loop over GLMs'
            'from different jobs and subjects.'}';

function vout = nirs_cfg_vout_AnalyzeGLM(job)
vout = cfg_dep;                     
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROC - Receiver Operating Curve Module - Sensitivity and specificity test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ROCDeleteLarge      = cfg_menu;
ROCDeleteLarge.tag  = 'ROCDeleteLarge';
ROCDeleteLarge.name = 'Delete large files';
ROCDeleteLarge.labels = {'Yes','No'};
ROCDeleteLarge.values = {1,2};
ROCDeleteLarge.val = {1};
ROCDeleteLarge.help = {'Delete large files (.nir) after each estimation.'};

ROCiternum         = cfg_entry; 
ROCiternum.name    = 'Number of iterations';
ROCiternum.tag     = 'ROCiternum';  
ROCiternum.val{1}  = 10;
ROCiternum.strtype = 'r';
ROCiternum.num     = [1 1];     
ROCiternum.help    = {'Number of iterations'}; 

% Executable Branch
ROCtest      = cfg_exbranch;       
ROCtest.name = 'ROC Sensitivity and specificity testing';            
ROCtest.tag  = 'ROCtest';
ROCtest.val  = {NIRSmat ROCLoopJob ROCDeleteLarge ROCiternum}; 
ROCtest.prog = @nirs_run_ROCtest;  
ROCtest.vout = @nirs_cfg_vout_ROCtest; 
ROCtest.help = {'This module performs a large loop over GLMs'
            'specified with different random seeds. '
            'To use it, user need to first specify in the Matlabbatch'
            'front end a sequence of modules to be run, typically starting'
            'with a module that requires a random seed (such as the '
            'AddTestStimuli module). This sequence of modules is referred '
            'to as a job. The code will run that job repetitively '
            'by incrementing the random seed as many times as specified. '}';

function vout = nirs_cfg_vout_ROCtest(job)
vout = cfg_dep;                     
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NIRS_SPM GLM Results Display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Executable Branch
NIRS_SPM_model_display      = cfg_exbranch;       
NIRS_SPM_model_display.name = 'NIRS_SPM Results Display';            
NIRS_SPM_model_display.tag  = 'NIRS_SPM_model_display';
NIRS_SPM_model_display.val  = {NIRS_SPM_Coregistration_Channels}; 
NIRS_SPM_model_display.prog = @nirs_run_NIRS_SPM_model_display;  
NIRS_SPM_model_display.vout = @nirs_cfg_vout_NIRS_SPM_model_display; 
NIRS_SPM_model_display.help = {'NIRS_SPM Results Display.'};

function vout = nirs_cfg_vout_NIRS_SPM_model_display(job)
vout = cfg_dep;                     
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NIRS_SPM Contrasts Display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Select map of interest
map_file         = cfg_files;  
map_file.name    = 'Select statistical map'; % The displayed name
map_file.tag     = 'map_file';       
map_file.filter  = 'mat';    
map_file.num     = [1 1];     % Number of inputs required 
map_file.help    = {'Select statistical map of interest (file '
            'containing SPM_nirs (for group) or cinterp_SPM_nirs '
            'structure (for individual) ) - used to select channel '
            'with nearest projected distance to statistical map maximum.'}'; % help text displayed

%select data file of interest
data_file         = cfg_files;  
data_file.name    = 'Select NIRS_SPM data file'; 
data_file.tag     = 'data_file';       
data_file.filter  = 'mat';    
data_file.num     = [1 1];     
data_file.help    = {'Select data file to calculate contrast estimates '
                '(file containing SPM_nirs structure) -- does not have '
                'to be related to statistical map.'}'; 

%Select contrasts
reg_num         = cfg_entry; %
reg_num.name    = 'Regressor identification numbers'; 
reg_num.tag     = 'reg_num';       
reg_num.strtype = 'r'; 
reg_num.num     = [1 Inf];     
reg_num.help    = {'Enter regressor numbers as a Matlab row vector '
        '(get from the design matrix associated with the data file)'}'; 

% Executable Branch
NIRS_SPM_contrast_display      = cfg_exbranch;       
NIRS_SPM_contrast_display.name = 'NIRS_SPM Contrasts Estimate Display';             
NIRS_SPM_contrast_display.tag  = 'NIRS_SPM_contrast_display'; 
%NIRSmat not used currently
NIRS_SPM_contrast_display.val  = {map_file data_file ...
            NIRS_SPM_Coregistration_Channels view reg_num};  
NIRS_SPM_contrast_display.prog = @nirs_run_NIRS_SPM_contrast_display;  
NIRS_SPM_contrast_display.vout = @nirs_cfg_vout_NIRS_SPM_contrast_display; 
NIRS_SPM_contrast_display.help = {'NIRS_SPM Contrast Estimates Results Display.'};

function vout = nirs_cfg_vout_NIRS_SPM_contrast_display(job)
vout = cfg_dep;                     
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration NIRS_SPM diagnostics for protocole and detrending time-series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%select data files of interest
data_files         = cfg_files;  
data_files.name    = 'Select NIRS_SPM data files'; 
data_files.tag     = 'data_files';       
data_files.filter  = 'mat';    
data_files.num     = [1 Inf];    
data_files.help    = {'Select data files for diagnostic (containing '
                'SPM_nirs structure).'}'; 

ch_num         = cfg_entry; 
ch_num.name    = 'Channel identification numbers'; 
ch_num.tag     = 'ch_num';       
ch_num.strtype = 'r';     
ch_num.num     = [1 Inf];     
ch_num.help    = {'Enter channel numbers as a Matlab row vector '
            '(get from the design matrix associated with the data file)'}'; 

% Executable Branch
NIRS_SPM_diagnostic      = cfg_exbranch;       
NIRS_SPM_diagnostic.name = 'NIRS_SPM diagnostics for protocole and detrending time-series';             
NIRS_SPM_diagnostic.tag  = 'NIRS_SPM_diagnostic'; 
NIRS_SPM_diagnostic.val  = {data_files ch_num}; 
NIRS_SPM_diagnostic.prog = @nirs_run_NIRS_SPM_diagnostic; 
NIRS_SPM_diagnostic.vout = @nirs_cfg_vout_NIRS_SPM_diagnostic; 
NIRS_SPM_diagnostic.help = {'NIRS_SPM diagnostics for protocole and detrending time-series.'};

function vout = nirs_cfg_vout_NIRS_SPM_diagnostic(job)
vout = cfg_dep;                     
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration NIRS Hemodynamic Modeling HDM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spmmat         = cfg_files; %Select NIRS.mat for this subject 
spmmat.name    = 'Select SPM.mat (for BOLD)'; % The displayed name
spmmat.tag     = 'spmmat';       %file names
spmmat.filter = 'mat';
spmmat.ufilter = '^SPM.mat$';    
spmmat.num     = [1 1];     % Number of inputs required 
spmmat.help    = {'Select SPM.mat of BOLD estimation for this subject.'}; 

spmmat_ASL         = cfg_files; 
spmmat_ASL.name    = 'Select SPM.mat (for ASL)'; % The displayed name
spmmat_ASL.tag     = 'spmmat_ASL';       %file names
spmmat_ASL.filter = 'mat';
spmmat_ASL.ufilter = '^SPM.mat$';    
spmmat_ASL.num     = [1 1];     % Number of inputs required 
spmmat_ASL.help    = {'Select SPM.mat of ASL estimation for this subject.'}; 

BOLD           = cfg_branch;
BOLD.name      = 'BOLD only';
BOLD.tag       = 'BOLD';
BOLD.val       = {spmmat}; 
BOLD.help      = {''};

ASL           = cfg_branch;
ASL.name      = 'ASL only';
ASL.tag       = 'ASL';
ASL.val       = {spmmat_ASL}; 
ASL.help      = {''};

BOLD_ASL           = cfg_branch;
BOLD_ASL.name      = 'BOLD and ASL';
BOLD_ASL.tag       = 'BOLD_ASL';
BOLD_ASL.val       = {spmmat spmmat_ASL}; 
BOLD_ASL.help      = {''};

spmmat_HbO         = cfg_files; 
spmmat_HbO.name    = 'Select SPM.mat (for HbO)'; % The displayed name
spmmat_HbO.tag     = 'spmmat_HbO';       %file names
spmmat_HbO.filter = 'mat';   
spmmat_HbO.num     = [1 1];     % Number of inputs required 
spmmat_HbO.help    = {'Select SPM_nirs matrix of HbO estimation for this subject.'};

spmmat_HbR         = cfg_files; 
spmmat_HbR.name    = 'Select SPM.mat (for HbR)'; % The displayed name
spmmat_HbR.tag     = 'spmmat_HbR';       %file names
spmmat_HbR.filter = 'mat';
spmmat_HbR.num     = [1 1];     
spmmat_HbR.help    = {'Select SPM_nirs matrix of HbR estimation for this subject.'}; 

ch_num_avg         = cfg_entry;
ch_num_avg.name    = 'Channel identification numbers';
ch_num_avg.tag     = 'ch_num_avg';       
ch_num_avg.strtype = 'r';
ch_num_avg.num     = [1 Inf];   
ch_num_avg.help    = {'Enter channel numbers to be averaged over as a '
                    'Matlab row vector, instead of map_file.'}'; 
                
stat_map_mode      = cfg_branch;
stat_map_mode.name      = 'Stat map channel selection mode';
stat_map_mode.tag       = 'stat_map_mode';
stat_map_mode.val       = {map_file}; 
stat_map_mode.help      = {'Channel selection via maximum of statistical map.'};

ch_avg_mode      = cfg_branch;
ch_avg_mode.name      = 'Channel selection by average';
ch_avg_mode.tag       = 'ch_avg_mode';
ch_avg_mode.val       = {ch_num_avg}; 
ch_avg_mode.help      = {'Channel selection via maximum of statistical map.'};

ch_mode           = cfg_choice;
ch_mode.name      = 'Channel Selection Mode';
ch_mode.tag       = 'ch_mode';
ch_mode.values    = {stat_map_mode ch_avg_mode}; 
ch_mode.val       = {ch_avg_mode}; 
ch_mode.help      = {'Choose channel selection mode: (1) take max of a '
        'statistical map (2) specify channels to be averaged over.'}'; 

downsamplingFactor1      = cfg_entry;
downsamplingFactor1.tag  = 'downsamplingFactor1';
downsamplingFactor1.name = 'Downsampling Factor';
downsamplingFactor1.val = {5};
downsamplingFactor1.strtype = 'r';  
downsamplingFactor1.num     = [1 1]; 
downsamplingFactor1.help    = {'Specify downsampling factor of time series.'}; 

rescalingFactor1      = cfg_entry;
rescalingFactor1.tag  = 'rescalingFactor1';
rescalingFactor1.name = 'Rescaling Factor';
rescalingFactor1.val = {0.05};
rescalingFactor1.strtype = 'r';  
rescalingFactor1.num     = [1 1]; 
rescalingFactor1.help    = {'Specify rescaling factor.'}; 

HbO_HbR            = cfg_branch;
HbO_HbR.name      = 'HbO+HbR';
HbO_HbR.tag       = 'HbO_HbR';
HbO_HbR.val       = {ch_mode spmmat_HbO spmmat_HbR NIRS_SPM_Coregistration_Channels view downsamplingFactor1 rescalingFactor1}; 
HbO_HbR.help      = {'HbO and HbR HDM estimation'};

Modalities           = cfg_choice;
Modalities.name      = 'Modalities: BOLD, BOLD + ASL estimation, ASL only, HbO+HbR:';
Modalities.tag       = 'Modalities';
%Modalities.labels    = {'BOLD' 'BOLD + ASL' 'ASL' 'HbO+HbR'};
Modalities.values    = {BOLD BOLD_ASL ASL HbO_HbR}; %{BOLDchoice spmmat_ASL}; %{'BOLD' 'BOLD + ASL'};
Modalities.val       = {HbO_HbR}; %{BOLD_ASL}; %{spmmat_ASL};
Modalities.help      = {'Choose data type: BOLD, BOLD + ASL, ASL only, HbO+HbR'}; %, HbO+HbR, HbO+HbR+speckle for hemodynamic parameters estimation'};

Model_Choice           = cfg_menu;
Model_Choice.name      = 'Choice of Model';
Model_Choice.tag       = 'Model_Choice';
Model_Choice.labels    = {'Buxton-Friston' 'Zheng-Mayhew' 'Huppert1'};
Model_Choice.values    = {0,1,2};
Model_Choice.val       = {0};
%Model_Choice.def  = @(val)nirs_get_defaults('readNIRS.boxy1.save_bin_dot', val{:});
Model_Choice.help      = {'Choose hemodynamic model: Buxton-Friston, '
    'Zheng-Mayhew, or 1-Compartment Huppert Model'}';

Stimuli     = cfg_entry; 
Stimuli.name    = 'Stimuli identification numbers'; 
Stimuli.tag     = 'Stimuli';       
Stimuli.strtype = 'r'; 
Stimuli.val     = {1};
Stimuli.num     = [1 Inf];     
Stimuli.help    = {'Enter stimuli numbers to include as a Matlab row '
    'vector (get from the design matrix associated with the data file)'}'; 


% Executable Branch
NIRS_HDM      = cfg_exbranch;       
NIRS_HDM.name = 'NIRS Hemodynamic Modelling';             
NIRS_HDM.tag  = 'NIRS_HDM'; 
NIRS_HDM.val  = {Modalities Model_Choice Stimuli}; 
NIRS_HDM.prog = @nirs_run_NIRS_HDM;  
NIRS_HDM.vout = @nirs_cfg_vout_NIRS_HDM; 
NIRS_HDM.help = {'NIRS_SPM Hemodynamic Modeling.'};

function vout = nirs_cfg_vout_NIRS_HDM(job)
vout = cfg_dep;                     
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Module 13 : CRIUGM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

heart_pace   = cfg_menu;
heart_pace.tag  = 'heart_pace';
heart_pace.name = 'Calculating heart pace ?';
heart_pace.labels = {'Yes','No'};
heart_pace.values = {1,0};
heart_pace.def  = @(val)nirs_get_defaults('readNIRS.criugm1.heart_pace', val{:});
heart_pace.help = {'If Yes some processings will be changed.'};

runVOIRE1      = cfg_exbranch;
runVOIRE1.name = 'Run VOIRE analysis';
runVOIRE1.tag  = 'runVOIRE1';
runVOIRE1.val  = {NIRSmat heart_pace criugm_paces1};
runVOIRE1.prog = @nirs_run_runVOIRE;
runVOIRE1.vout = @nirs_cfg_vout_runVOIRE;
runVOIRE1.help = {'.'};

    function vout = nirs_cfg_vout_runVOIRE(job)
        vout = cfg_dep;
        vout.sname      = 'NIRS.mat';
        vout.src_output = substruct('.','NIRSmat');
        vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Module 14 : CRIUGM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

acc_file         = cfg_files;
acc_file.name    = 'Accelerometer file'; % The displayed name
acc_file.tag     = 'acc_file';       %file names
acc_file.filter  = 'csv';
acc_file.num     = [1 Inf];     % Number of inputs required
acc_file.help    = {''}; % help text displayed

subj         = cfg_branch;
subj.tag     = 'subj';
subj.name    = 'Subject';
subj.val     = {acc_file};
subj.help    = {'Subject'};

generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Subjects';
generic.help    = {'Help'};
generic.values  = {subj};
generic.num     = [1 Inf];

runMOB1      = cfg_exbranch;
runMOB1.name = 'Run MOB analysis';
runMOB1.tag  = 'runMOB1';
runMOB1.val  = {generic};
runMOB1.prog = @nirs_run_runMOB;
runMOB1.vout = @nirs_cfg_vout_runMOB;
runMOB1.help = {'.'};

    function vout = nirs_cfg_vout_runMOB(job)
        vout = cfg_dep;
        vout.sname      = 'NIRS.mat';
        vout.src_output = substruct('.','NIRSmat');
        vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration main modules  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
%module 1
readNIRS        = cfg_choice; 
readNIRS.name   = 'Read NIRS data';
readNIRS.tag    = 'readNIRS';
readNIRS.values = {boxy1 criugm1 lot1}; 
readNIRS.help   = {'These modules read NIRS data in different formats.'};

%module 0: utilities to read onsets and create GLM stimuli structure
readOnsets        = cfg_choice; 
readOnsets.name   = 'Read Onsets';
readOnsets.tag    = 'readOnsets';
readOnsets.values = {AnalyzerOnsets readEprimeOnsets permuteOnsets addTestStimuli}; 
readOnsets.help   = {'These modules create stimuli structures '
                    'as inputs to the General Linear Model.'}';
%module 2
preprocANAT        = cfg_choice;
preprocANAT.name   = 'Preprocess anatomical image';
preprocANAT.tag    = 'preprocANAT';
preprocANAT.values = {detectVitamins1 MCsegment1 buildroi1}; 
preprocANAT.help   = {'These modules pre-process anatomical images '
        'so that clean anatomical images can be used with functional data.'}';

%module 3
coregNIRS        = cfg_choice; %cfg_repeat; 
coregNIRS.name   = 'Coregister NIRS data';
coregNIRS.tag    = 'coregNIRS';
coregNIRS.values = {coreg1 coreg_manual1 view3d1 resize1};
coregNIRS.help   = {'These modules perform coregistration ',...
            'between NIRS and an anatomical image.'};

%module 4 - NIRS preprocessing (heart rate detection, pruning bad channels, filters)
preprocessNIRS        = cfg_choice; 
preprocessNIRS.name   = 'Preprocess NIRS data';
preprocessNIRS.tag    = 'preprocessNIRS';
preprocessNIRS.values = {remove_chn_stdev criugm_paces1  ...
         mark_movement normalize_baseline ODtoHbOHbR generate_vhdr_vmrk}; %mark_negative HPF_LPF
preprocessNIRS.help   = {'These modules preprocess NIRS data: '
    'heart rate check, '
    'downsampling, removal of bad channels, filters.'}';

%module 8 - reconstruction
model_reconstruct    = cfg_choice; %cfg_repeat;
model_reconstruct.name   = '3D Reconstruction of NIRS data';
model_reconstruct.tag    = 'model_reconstruct';
model_reconstruct.values = {tikhonov1 ReMLreconstruct1};
model_reconstruct.help   = {'3D Reconstruction of NIRS data.'};

%module 9
model_specify        = cfg_choice; %cfg_repeat;
model_specify.name   = 'GLM Specification';
model_specify.tag    = 'model_specify';
model_specify.values = {wls_bglm_specify}; %NIRS_SPM_specify NIRS_SPM_specify_batch
model_specify.help   = {'These modules specify a GLM.'};

%module 10
model_estimate        = cfg_choice; %cfg_repeat; 
model_estimate.name   = 'GLM Estimation';
model_estimate.tag    = 'model_estimate';
model_estimate.values = {wls_bglm_estimate liom_contrast  ...
            liom_group AnalyzeGLM ROCtest}; 
        %NIRS_SPM_HPF_LPF NIRS_SPM_estimate NIRS_SPM_estimate_batch
        %NIRS_SPM_contrast NIRS_SPM_group
model_estimate.help   = {'These modules estimate a GLM.'};
 
%module 11
model_display        = cfg_choice; %cfg_repeat; 
model_display.name   = 'GLM Results Display';
model_display.tag    = 'model_display';
model_display.values = {NIRS_SPM_model_display ...
        NIRS_SPM_contrast_display NIRS_SPM_diagnostic};
model_display.help   = {'Display results of NIRS GLM.'};

%module 13
CRIUGM        = cfg_choice; %cfg_repeat; 
CRIUGM.name   = 'CRIUGM';
CRIUGM.tag    = 'CRIUGM';
CRIUGM.values = {runVOIRE1 runMOB1};
CRIUGM.help   = {'Data analysis for CRIUGM projects'};

%-----------------------------------------------------------------------
nirs10        = cfg_choice;
nirs10.name   = 'nirs10';
nirs10.tag    = 'nirs10'; %Careful, this tag nirs10 must be the same as
%the name of the toolbox and when called by spm_jobman in nirs10.m
nirs10.values = {readNIRS readOnsets preprocessNIRS preprocANAT coregNIRS ...
    configMC1 runMC1 makesens1 model_reconstruct model_specify ...
    model_estimate NIRS_HDM CRIUGM}; %model_display
end