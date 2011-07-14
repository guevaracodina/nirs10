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
age1.def     = @(val)nirs_get_defaults('readNIRS.criugm1.generic2.subj.age1', val{:});
age1.help    = {'Age of the subject. Used later for OD to HbO/HbR conversion.'};

text_brainsight         = cfg_files;
text_brainsight.tag     = 'text_brainsight';
text_brainsight.name    = 'Text file from Brainsight';
text_brainsight.filter  = '.txt';
text_brainsight.ufilter = '.*';
text_brainsight.num     = [1 1];
text_brainsight.help    = {'Select the text file from Brainsight.'};
 
T1_vitamins      = cfg_branch; 
T1_vitamins.name = 'Vitamins markers on T1';
T1_vitamins.tag  = 'T1_vitamins';
T1_vitamins.help = {'The helmet will be read in future module from T1 image and positions of vitamins in the image.'};

no_helmet      = cfg_branch; 
no_helmet.name = 'No helmet information';
no_helmet.tag  = 'no_helmet';
no_helmet.help = {'Helmet informations will be extracted from ''.nirs'' file.'};

helm_temp         = cfg_files;
helm_temp.tag     = 'helm_temp';
helm_temp.name    = 'Helmet template';
helm_temp.help = {'If you have chosen before template in choice : ''Individual T1 or template''.'};
helm_temp.filter  = 'dir';
helm_temp.ufilter = '.*';
helm_temp.val{1} = {''};
helm_temp.num     = [0 0];

% helmet         = cfg_choice;
% helmet.tag     = 'helmet';
% helmet.name    = 'Helmet';
% if indvdata.val=={indvdata_chosen}
%     helmet.values = {helm_temp};
% else
%     helmet.values = {text_brainsight T1_vitamins no_helmet};
% end
% helmet.val     = {helm_temp};
% helmet.help    = {'If you choose a Brainsight text file, it will be used to determine all you need about sources, detectors and other points of interest.'};

nirs_files         = cfg_files;
nirs_files.name    = '''^.nirs'' files';
nirs_files.tag     = 'nirs_files';
nirs_files.filter  = 'nirs';
nirs_files.num     = [0 Inf];
nirs_files.val{1} = {''};
nirs_files.help    = {'Select all the sessions sharing the same device and helmet. If no ''^.nirs'' is selected you won''t be able to run coregistration.'};

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

boldmask        = cfg_files;
boldmask.tag    = 'boldmask';
boldmask.name   = 'BOLD mask';
boldmask.filter  = 'image';  
boldmask.ufilter = '.*';
boldmask.num    = [0 Inf];
boldmask.val{1} = {''};
boldmask.help   = {'Help'};

subj_id         = cfg_entry;
subj_id.tag     = 'subj_id';
subj_id.name    = 'Subject ID';
subj_id.strtype = 's';
subj_id.num     = [1 Inf];
subj_id.val     = {}; 
subj_id.help    = {'A number must be entered.'};

choose_path         = cfg_entry;
choose_path.tag     = 'choose_path';
choose_path.name    = 'Path for study being created';
choose_path.strtype = 's';
choose_path.num     = [1 Inf];
choose_path.val     = {}; 
choose_path.help    = {'.'};

existing_study        = cfg_files;
existing_study.tag     = 'existing_study';
existing_study.name    = 'Existing study';
existing_study.filter  = 'dir';
existing_study.ufilter = '.*';
existing_study.num     = [1 1];
existing_study.help    = {'Choose directory where you want to put your study. If the directory does not exist you must create it then choose it in the batch.'};

study_path           = cfg_choice;
study_path.name      = 'Study path configuration';
study_path.tag       = 'study_path';
study_path.values    = {existing_study choose_path};%choose_path
study_path.val       = {existing_study}; 
study_path.help      = {'Choose the study the subject belongs to or specify a path (entire name should look like .\study_name).'}; 

indvdata_chosen      = cfg_branch; 
indvdata_chosen.name = 'One set of data per subject';
indvdata_chosen.tag  = 'indvdata_chosen';
indvdata_chosen.help = {'You will have to choose one T1 image in the field ''Raw anatomical image'' and one helmet in the field ''Helmet->Text file from Brainsight.''.'};

anatT1_template         = cfg_files; 
anatT1_template.name    = 'Anatomical template image'; 
anatT1_template.tag     = 'anatT1_template';       %file names
anatT1_template.filter  = 'image';  
anatT1_template.ufilter = '.*';
anatT1_template.def = @(val)nirs_get_defaults('coregNIRS.coreg1.anatT1_template', val{:});
anatT1_template.num     = [1 1];     % Number of inputs required 
anatT1_template.help    = {'Select anatomical template image for this subject.'};

anatT1_subj0         = cfg_files; 
anatT1_subj0.name    = 'Anatomical image of subject 0'; 
anatT1_subj0.tag     = 'anatT1_subj0';
anatT1_subj0.filter  = 'image';  
anatT1_subj0.ufilter = '.*';
anatT1_subj0.num     = [1 1];
anatT1_subj0.help    = {'Select anatomical template image for the subject 0.'};

template_chosen      = cfg_branch;
template_chosen.tag  = 'template_chosen';
template_chosen.name = 'Template';
template_chosen.val  = {anatT1_subj0 text_brainsight}; 
template_chosen.help = {'You must have a T1 image and a Brainsight registration of the right helmet for one subject (called subject 0).'};

indvdata        = cfg_choice;
indvdata.name   = 'Individual data or template for all';
indvdata.tag    = 'indvdata';
indvdata.values = {template_chosen indvdata_chosen};%choose_path
indvdata.val    = {indvdata_chosen};
indvdata.help   = {['Individual data allows you to choose data for each of the subject.'...
    'Template for all allows you to use one coregistration for all your subjects. You will need one T1 image and the registration of the helmet on the same person.']}; 



helmet         = cfg_choice;
helmet.tag     = 'helmet';
helmet.name    = 'Helmet';
helmet.values  = {text_brainsight T1_vitamins no_helmet helm_temp};
helmet.val     = {helm_temp};
helmet.help    = {'If you choose a Brainsight text file, it will be used to determine all you need about sources, detectors and other points of interest.'};

subj         = cfg_branch;
subj.tag     = 'subj';
subj.name    = 'Subject';
subj.val     = {subj_id age1 anatT1 helmet CWsystem nirs_files protocol TopoData boldmask};%config_path2
subj.help    = {'Subject'};

generic2         = cfg_repeat;
generic2.tag     = 'generic2';
generic2.name    = 'Subjects';
generic2.help    = {'Help'};
generic2.values  = {subj};
generic2.num     = [1 Inf];

study_cfg         = cfg_branch;
study_cfg.tag     = 'study_cfg';
study_cfg.name    = 'Study configuration';
study_cfg.val     = {study_path indvdata}; 
study_cfg.help    = {''};

% The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
% Executable Branch
criugm1      = cfg_exbranch;
criugm1.name = 'Read and format CRIUGM data';
criugm1.tag  = 'criugm1';
criugm1.val  = {study_cfg  generic2};
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

testGamma         = cfg_menu;
testGamma.tag     = 'testGamma';
testGamma.name    = 'Use Gamma function';
testGamma.help    = {
                'Use gamma function or canonical HRF for test'
                'Note that canonical HRF is leads to lower estimability due to the undershoot'}';
testGamma.labels = {
               'Gamma Function'
               'Canonical HRF'
                }';
testGamma.values = {1, 2};
testGamma.val = {1};

testFilterX         = cfg_menu;
testFilterX.tag     = 'testFilterX';
testFilterX.name    = 'Filter design matrix';
testFilterX.help    = {
                'Filter design matrix prior to calculating power of protocol'}';
testFilterX.labels = {
               'Filter On'
               'Filter Off'
                }';
testFilterX.values = {1, 0};
testFilterX.val = {1};

testFilterData         = cfg_menu;
testFilterData.tag     = 'testFilterData';
testFilterData.name    = 'Filter data';
testFilterData.help    = {
                'Filter data channels prior to calculating power of baseline'}';
testFilterData.labels = {
               'Filter On'
               'Filter Off'
                }';
testFilterData.values = {1, 0};
testFilterData.val = {1};

testStdvsPower        = cfg_menu;
testStdvsPower.tag     = 'testStdvsPower';
testStdvsPower.name    = 'Normalization choice';
testStdvsPower.help    = {
                'Normalize added protocol based on channel standard deviation or power'
                'Either option should give the same result...'}';
testStdvsPower.labels = {
               'Standard deviation'
               '(square root of) Power'
                }';
testStdvsPower.values = {1, 0};
testStdvsPower.val = {1};


testBfNorm        = cfg_menu;
testBfNorm.tag     = 'testBfNorm';
testBfNorm.name    = 'Further Normalization choice';
testBfNorm.help    = {
                'Normalize further by rescaling both 1st and 2nd Volterra by the max value of the 1st Volterra'}';
testBfNorm.labels = {
               'Normalize'
               'Do not Normalize'
                }';
testBfNorm.values = {1, 0};
testBfNorm.val = {0};

testHPFButterOn         = cfg_menu;
testHPFButterOn.tag     = 'testHPFButterOn';
testHPFButterOn.name    = 'Butter HPF';
testHPFButterOn.help    = {
                'Preferred option: ON, at 0.004'}';
testHPFButterOn.labels = {
               'Filter On'
               'Filter Off'
                }';
testHPFButterOn.values = {1, 0};
testHPFButterOn.val = {1};

testHPFbutterCutoff         = cfg_entry; 
testHPFbutterCutoff.name    = 'HPF Butter cutoff';
testHPFbutterCutoff.tag     = 'testHPFbutterCutoff';       
testHPFbutterCutoff.strtype = 'r';
testHPFbutterCutoff.num     = [1 1];     
testHPFbutterCutoff.val = {0.004};
testHPFbutterCutoff.help    = {'HPF cutoff in Hz. (0.004 Hz recommended)'}';

testHPFbutterOrder         = cfg_entry; 
testHPFbutterOrder.name    = 'HPF Butter order';
testHPFbutterOrder.tag     = 'testHPFbutterOrder';       
testHPFbutterOrder.strtype = 'r';
testHPFbutterOrder.num     = [1 1];     
testHPFbutterOrder.val = {5};
testHPFbutterOrder.help    = {'HPF order (3 recommended)'}';

testLPFGaussianOn         = cfg_menu;
testLPFGaussianOn.tag     = 'testLPFGaussianOn';
testLPFGaussianOn.name    = 'Gaussian LPF';
testLPFGaussianOn.help    = {
                'Preferred option: ON, at 1.5 s'}';
testLPFGaussianOn.labels = {
               'Filter On'
               'Filter Off'
                }';
testLPFGaussianOn.values = {1, 0};
testLPFGaussianOn.val = {1};

testLPFGaussianFWHM         = cfg_entry; 
testLPFGaussianFWHM.name    = 'LPF Gaussian FWHM';
testLPFGaussianFWHM.tag     = 'testLPFGaussianFWHM';       
testLPFGaussianFWHM.strtype = 'r';
testLPFGaussianFWHM.num     = [1 1];     
testLPFGaussianFWHM.val = {1.5};
testLPFGaussianFWHM.help    = {'LPF Gaussian FWHM in seconds ( 1.5 s recommended)'}';

testWaveletMDLOn         = cfg_menu;
testWaveletMDLOn.tag     = 'testWaveletMDLOn';
testWaveletMDLOn.name    = 'Wavelet MDL HPF';
testWaveletMDLOn.help    = {
                'Preferred option: Off'}';
testWaveletMDLOn.labels = {
               'Filter On'
               'Filter Off'
                }';
testWaveletMDLOn.values = {1, 0};
testWaveletMDLOn.val = {0};


% Executable Branch
addTestStimuli      = cfg_exbranch;       
addTestStimuli.name = 'Add Stimuli with HRFs for testing';             
addTestStimuli.tag  = 'addTestStimuli'; 
addTestStimuli.val  = {NIRSmat DelPreviousData NewDirCopyNIRS testStimulusName testStimuliNumber ...
                testSessionNumber testWavelength testAmplitudeTarget ...
                voltAddStim testAmplitude2 keepAllChannels testChannels testDupChannels testPType ...
                testGamma testFilterX testFilterData testStdvsPower testBfNorm ...
                testHPFButterOn testHPFbutterCutoff testHPFbutterOrder ...
                testLPFGaussianOn testLPFGaussianFWHM testWaveletMDLOn};   
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
image_in.val{1} = {''};
image_in.num     = [0 1];

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
buildroi1.name = 'Build ROI';
buildroi1.val  = {NIRSmat NewDirCopyNIRS keepAllChannels image_in output_prefix};
buildroi1.prog = @nirs_run_buildroi2;
buildroi1.vout = @nirs_cfg_vout_buildroi1;
buildroi1.help = {'Define region of interest containing all the selected channels. Please only enter the channels numbers for the first wavelength.'};

function vout = nirs_cfg_vout_buildroi1(job)
vout = cfg_dep;                     
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

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
    rebel_surrounding rebel_thresh_hs process_image}; %%%%%% pk NIRSmat_optional ???
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

View6Projections      = cfg_menu;
View6Projections.tag  = 'View6Projections';
View6Projections.name = 'View the 6 Projections';
View6Projections.labels = {'True','False'};
View6Projections.values = {1,0};
View6Projections.val  = {0};
View6Projections.help = {'View channel locations for the 6 projections.'}';

render_file         = cfg_files; 
render_file.name    = 'Render file'; 
render_file.tag     = 'render_file';      
%render_file.filter  = 'image';  
%render_file.ufilter = '.*';
render_file.num     = [1 1];     
render_file.help    = {'Grey matter (c1) anatomical image or rendered version of this c1 image.'
    'Normalized or not according to Normalization choice option below.'}';
        
        
render_normalize_choice = cfg_menu;
render_normalize_choice.tag    = 'render_normalize_choice';
render_normalize_choice.name   = 'Normalization Choice';
render_normalize_choice.labels = {'MNI Talairach Tournoux atlas','Subject Coordinates'};
render_normalize_choice.values = {1,0};
render_normalize_choice.val = {1};
render_normalize_choice.help   = {'Normalization choice. Need to specify appropriate'
    'file: normalized file if normalized to TT, unnormalized if in subject coordinates.'
    'Method in subject coordinates not coded up yet.'}';
        
render_template         = cfg_branch;
render_template.tag     = 'render_template';
render_template.name    = 'Render to SPM single subject template';
render_template.val     = {}; 
render_template.help    = {'Render to template.'};

render_subject         = cfg_branch;
render_subject.tag     = 'render_subject';
render_subject.name    = 'Render to subject';
render_subject.val     = {render_file render_normalize_choice}; 
render_subject.help    = {'Render to subject. OPTION NOT FUNCTIONAL YET: '
    'Problem with coordinate systems and projections.'}';

render_choice        = cfg_choice;
render_choice.name   = 'Render to template or subject';
render_choice.tag    = 'render_choice';
render_choice.values = {render_template,render_subject};
render_choice.val    = {render_template};
render_choice.help   = {'Render to template or subject.'};

coreg1      = cfg_exbranch;       
coreg1.name = 'NIRScoreg';             
coreg1.tag  = 'coreg1'; 
coreg1.val  = {NIRSmat DelPreviousData NewDirCopyNIRS anatT1 segT1_4fit ...
    anatT1_template fid_in_subject_MNI nasion_wMNI AL_wMNI AR_wMNI ...
    GenDataTopo render_choice View6Projections};    
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
% Coreg vers le template T1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coreg2      = cfg_exbranch;       
coreg2.name = 'NIRScoreg with T1 template';             
coreg2.tag  = 'coreg2'; 
coreg2.val  = {NIRSmat DelPreviousData NewDirCopyNIRS anatT1 segT1_4fit ...
    anatT1_template fid_in_subject_MNI nasion_wMNI AL_wMNI AR_wMNI GenDataTopo};    
coreg2.prog = @nirs_run_coreg_2templateT1;  
coreg2.vout = @nirs_cfg_vout_coreg2; 
coreg2.help = {'Automatic coregistration with T1 template. Use this choice in the case you don''t have the anatomical T1 images of the subject.'};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_coreg2(job)
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
    'preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_exercise.STFT_param2.win_type2', val{:});
win_type2.help = {'Only Hanning for now.'};

win_width2         = cfg_entry;
win_width2.name    = 'Window width';
win_width2.tag     = 'win_width2';       
win_width2.strtype = 'r';
win_width2.num     = [1 1];
win_width2.def  = @(val)nirs_get_defaults(...
    'preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_exercise.STFT_param2.win_width2', val{:});
win_width2.help    = {'Window width in SECONDS'};

Nprobe2         = cfg_entry;
Nprobe2.name    = 'Number of probes';
Nprobe2.tag     = 'Nprobe2';       
Nprobe2.strtype = 'r';
Nprobe2.num     = [1 1];
Nprobe2.def  = @(val)nirs_get_defaults(...
    'preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_exercise.STFT_param2.Nprobe2', val{:});
Nprobe2.help    = {'Number of probes taken along the signal (power of 2)'};

fft_size2         = cfg_entry;
fft_size2.name    = 'FFT size';
fft_size2.tag     = 'fft_size2';       
fft_size2.strtype = 'r';
fft_size2.num     = [1 1];
fft_size2.def  = @(val)nirs_get_defaults(...
    'preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_exercise.STFT_param2.fft_size2', val{:});
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
%     'preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_exercise.detect_wavelength2', val{:}); 
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
%     'preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_exercise.MinHeartRate2', val{:}); 
% MinHeartRate2.help    = {'Enter minimum heart rate allowed for final detection in Hz.'}';
% 
% MaxHeartRate2         = cfg_entry; 
% MaxHeartRate2.name    = 'Maximum Heart Rate for Detection';
% MaxHeartRate2.tag     = 'MaxHeartRate2';       
% MaxHeartRate2.strtype = 'r';
% MaxHeartRate2.num     = [1 1];     
% MaxHeartRate2.def     = @(val)nirs_get_defaults(...
%     'preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_exercise.MaxHeartRate2', val{:}); 
% MaxHeartRate2.help    = {'Enter maximum heart rate allowed for final detection in Hz.'}';
% 
% InternalMinHeartRate2         = cfg_entry; 
% InternalMinHeartRate2.name    = 'Internal Minimum Heart Rate for Detection';
% InternalMinHeartRate2.tag     = 'InternalMinHeartRate2';       
% InternalMinHeartRate2.strtype = 'r';
% InternalMinHeartRate2.num     = [1 1];     
% InternalMinHeartRate2.def     = @(val)nirs_get_defaults(...
%     'preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_exercise.InternalMinHeartRate2', val{:}); 
% InternalMinHeartRate2.help    = {'Enter minimum heart rate allowed for detection in Hz, '
%                         '"internal", i.e. for check over small data windows for FFT.'}';
% 
% InternalMaxHeartRate2         = cfg_entry; 
% InternalMaxHeartRate2.name    = 'Internal Maximum Heart Rate for Detection';
% InternalMaxHeartRate2.tag     = 'InternalMaxHeartRate2';       
% InternalMaxHeartRate2.strtype = 'r';
% InternalMaxHeartRate2.num     = [1 1];     
% InternalMaxHeartRate2.def     = @(val)nirs_get_defaults(...
%     'preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_exercise.InternalMaxHeartRate2', val{:}); 
% InternalMaxHeartRate2.help    = {'Enter maximum heart rate allowed for detection in Hz, '
%                         '"internal", i.e. for check over small data windows for FFT.'}';
% 
% MaxHeartStdev2         = cfg_entry; 
% MaxHeartStdev2.name    = 'Maximum Heart Rate Standard Deviation';
% MaxHeartStdev2.tag     = 'MaxHeartStdev2';       
% MaxHeartStdev2.strtype = 'r';
% MaxHeartStdev2.num     = [1 1];     
% MaxHeartStdev2.def     = @(val)nirs_get_defaults(...
%     'preprocessNIRS.criugm_paces1.heart_rate_cfg.heart_exercise.MaxHeartStdev2', val{:}); 
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
save_heart_rate_figure      = cfg_menu;
save_heart_rate_figure.tag  = 'save_heart_rate_figure';
save_heart_rate_figure.name = 'save_heart_rate_figure';
save_heart_rate_figure.labels = {'True','False'};
save_heart_rate_figure.values = {1,0};
save_heart_rate_figure.val  = {1};
save_heart_rate_figure.help = {'save_heart_rate_figure.'}';

display_heart_rate_figure      = cfg_menu;
display_heart_rate_figure.tag  = 'display_heart_rate_figure';
display_heart_rate_figure.name = 'display_heart_rate_figure';
display_heart_rate_figure.labels = {'True','False'};
display_heart_rate_figure.values = {1,0};
display_heart_rate_figure.val  = {0};
display_heart_rate_figure.help = {'display_heart_rate_figure.'}';

% Executable Branch
criugm_paces1      = cfg_exbranch;       
criugm_paces1.name = 'Heart rate utility';             
criugm_paces1.tag  = 'criugm_paces1'; 
criugm_paces1.val  = {NIRSmat DelPreviousData NewDirCopyNIRS heart_rate_cfg ...
    remove_no_heartbeat save_heart_rate_figure display_heart_rate_figure};   
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
nirs_hpf.val       = {hpf_none}; 
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


% % Executable Branch
% HPF_LPF      = cfg_exbranch;       
% HPF_LPF.name = 'Filters';             
% HPF_LPF.tag  = 'HPF_LPF'; 
% HPF_LPF.val  = {NIRSmat DelPreviousData NewDirCopyNIRS nirs_hpf nirs_lpf}; 
% HPF_LPF.prog = @nirs_run_HPF_LPF;  
% HPF_LPF.vout = @nirs_cfg_vout_HPF_LPF; 
% HPF_LPF.help = {'Filters: currently only low pass, with or without ',...
%     'a downsampling factor.'};
% 
% function vout = nirs_cfg_vout_HPF_LPF(job)
% vout = cfg_dep;                     
% vout.sname      = 'NIRS.mat';       
% vout.src_output = substruct('.','NIRSmat'); 
% vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
% end

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

mcim_in         = cfg_files;
mcim_in.name    = 'MC segmented volume';
mcim_in.tag     = 'mcim_in';
mcim_in.filter = 'image';
mcim_in.ufilter = '.nii';    
mcim_in.num     = [1 1];
mcim_in.help    = {'Select MC segmented volume for this subject.'};

mcim_cfg           = cfg_choice;
mcim_cfg.name      = 'Image';
mcim_cfg.tag       = 'mcim_cfg';
mcim_cfg.values    = {latest_mcim mcim_in};
mcim_cfg.val       = {latest_mcim}; 
mcim_cfg.help      = {'Choose latest ROI simulated or specify the segmented image.'}; 

MC_CUDAchoice    = cfg_menu;
MC_CUDAchoice.name   = 'Configuration file type';
MC_CUDAchoice.tag    = 'MC_CUDAchoice';
MC_CUDAchoice.labels = {'MCX: Monte Carlo Extreme','tMCimg','Both'};
MC_CUDAchoice.values = {1,2,3};
MC_CUDAchoice.def    = @(val)nirs_get_defaults('configMC1.MC_CUDAchoice', val{:});
MC_CUDAchoice.help   = {'Choose type of configuration files to generate.'};

no_pve      = cfg_branch; 
no_pve.name = 'No';
no_pve.tag  = 'no_pve';
no_pve.help = {'PVE won''t be calculated.'};

calc_pve         = cfg_files;
calc_pve.name    = 'Mask for PVE';
calc_pve.tag     = 'calc_pve';
calc_pve.filter  = 'image';
calc_pve.ufilter = '.nii';    
calc_pve.num     = [1 1];
calc_pve.help    = {'BOLD, ASL or any anatomical mask (from create mask module).'};

pve_cfg           = cfg_choice;
pve_cfg.name      = 'PVE';
pve_cfg.tag       = 'pve_cfg';
pve_cfg.values    = {no_pve calc_pve};
pve_cfg.val       = {calc_pve}; 
pve_cfg.help      = {'.'};

% est ce qu'il y a pas un pb du au fait qu'il attend un directory ???
MC_configdir         = cfg_entry;
MC_configdir.tag     = 'MC_configdir';
MC_configdir.name    = 'Monte Carlo configuration files directory';
MC_configdir.strtype = 's';
MC_configdir.num     = [1 Inf];
MC_configdir.def     = @(val)nirs_get_defaults('configMC1.MC_configdir', val{:}); 
MC_configdir.help    = {'Directory to put Monte Carlo configuration files.'
    'NO Longer USED'}';

MC_nam         = cfg_entry;
MC_nam.tag     = 'MC_nam';
MC_nam.name    = 'Monte Carlo simulation name';
MC_nam.strtype = 's';
MC_nam.num     = [1 Inf];
MC_nam.val{1}  = 'sim';
MC_nam.help    = {'Name of Monte Carlo simulation.'
    'If a simulation has already been run, '
    'the current date will be automatically added to the name specified'}';

%--------------------------------------------------------------------------
nphotons         = cfg_entry; 
nphotons.name    = 'Number of photons'; % The displayed name
nphotons.tag     = 'nphotons';       %file names
nphotons.strtype = 'r';  
nphotons.num     = [1 1];     % Number of inputs required 
nphotons.def = @(val)nirs_get_defaults('configMC1.nphotons', val{:});
nphotons.help    = {'Input number of photons (not currently used).'}; 
             
seed         = cfg_entry; 
seed.name    = 'Random seed';
seed.tag     = 'seed';
seed.strtype = 'r';  
seed.num     = [1 1];
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
configMC1.val  = {NIRSmat MC_nam mcim_cfg MC_CUDAchoice pve_cfg MC_configdir MC_parameters};    
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
MCXconfigFiles.num     = [0 Inf];      
MCXconfigFiles.val{1}  = {''};
MCXconfigFiles.help    = {'Select input files (.inp for MCX).'}; 

%%%%% Options for MC run
MCX_t      = cfg_entry;
MCX_t.tag  = 'MCX_t';
MCX_t.name = 'Thread number';
MCX_t.val = {4800};
MCX_t.strtype = 'r';  
MCX_t.num     = [1 1]; 
MCX_t.help = {'Total number of threads -- see the specifications of your GPU'
    'Examples: NVidia GeForce 570: 4800'
    'GeForce GTX295: 1792?'}';

MCX_T      = cfg_entry;
MCX_T.tag  = 'MCX_T';
MCX_T.name = 'Thread number per block';
MCX_T.val = {480};
MCX_T.strtype = 'r';  
MCX_T.num     = [1 1]; 
MCX_T.help = {'Blocksize -- see the specifications of your GPU'
    'Examples: NVidia GeForce 570: 480'
    'GeForce GTX295: ?'}';

MCX_r      = cfg_entry;
MCX_r.tag  = 'MCX_r';
MCX_r.name = 'Number of repetitions';
MCX_r.val = {1};
MCX_r.strtype = 'r';  
MCX_r.num     = [1 1]; 
MCX_r.help = {'Number of repetitions: number of times that the simulation'
    'will be repeated, with different random seeds, to increase the total'
    'number of photons'}';

MCX_g      = cfg_entry;
MCX_g.tag  = 'MCX_g';
MCX_g.name = 'Number of gates';
MCX_g.val = {1};
MCX_g.strtype = 'r';  
MCX_g.num     = [1 1]; 
MCX_g.help = {'Number of gates: if larger than the number of specified gates'
    'in the config files, then only the number of gates specified in the config'
    'files will be run. Otherwise, this allows the user to run with fewer gates.'}';

MCX_l      = cfg_menu;
MCX_l.tag  = 'MCX_l';
MCX_l.name = 'Write log file';
MCX_l.labels = {'Yes', 'No'};
MCX_l.values = {1,0};
MCX_l.val  = {1};
MCX_l.help = {'Write log file.'}'; 

MCX1         = cfg_branch;
MCX1.tag     = 'MCX1';
MCX1.name    = 'Monte Carlo Extreme';
MCX1.val     = {MCXconfigFiles MCX_t MCX_T MCX_r MCX_g MCX_l}; % MCXconfig};
MCX1.help    = {'Run Monte Carlo Extreme simulation'};

tMCimg_configFiles         = cfg_files; %Select 
tMCimg_configFiles.name    = 'Select input files'; 
tMCimg_configFiles.tag     = 'tMCimg_configFiles';       
tMCimg_configFiles.ufilter = '.cfg';    
tMCimg_configFiles.num     = [0 Inf];   
tMCimg_configFiles.val{1}  = {''};
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

MCtestOneChannel      = cfg_menu;
MCtestOneChannel.tag  = 'MCtestOneChannel';
MCtestOneChannel.name = 'Test MC by running only first channel';
MCtestOneChannel.labels = {'Yes','No'};
MCtestOneChannel.values = {1,0};
MCtestOneChannel.val  = {0};
MCtestOneChannel.help = {'To do a quick test, run simulation only on first source and first detector for one wavelength.'}';


% Executable Branch
runMC1      = cfg_exbranch;      
runMC1.name = 'Run Monte Carlo simulation';           
runMC1.tag  = 'runMC1';
runMC1.val  = {NIRSmat NewDirCopyNIRS MC_runCUDAchoice MCtestOneChannel}; 
runMC1.prog = @nirs_run_runMC;  
runMC1.vout = @nirs_cfg_vout_runMC; 
runMC1.help = {'Run Monte Carlo simulation.'};

%make .mc2 or (.his, .2pt) file names available as a dependency
function vout = nirs_cfg_vout_runMC(job)
vout = cfg_dep;                     
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
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

outMCfiles         = cfg_files;
outMCfiles.name    = 'Select MC output files';
outMCfiles.tag     = 'outMCfiles';
outMCfiles.ufilter = {'.2pt','.mc2'};    
outMCfiles.num     = [0 Inf];
outMCfiles.val     = {};
outMCfiles.help    = {'Select .mc2 or .2pt files for this subject.'}; 

% Executable Branch
makesens1      = cfg_exbranch;      
makesens1.name = 'Sensitivity Matrix';            
makesens1.tag  = 'makesens1'; 
makesens1.val  = {NIRSmat outMCfiles};
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
%4.8 Calculate partial volume effect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_in         = cfg_files;
dir_in.tag     = 'dir_in';
dir_in.name    = 'MonteCarlo output directory';
dir_in.help    = {'Select the MonteCarlo simulation output directory.'};
dir_in.filter = 'dir';
dir_in.val{1} = {''};
dir_in.ufilter = '.*';
dir_in.num     = [0 1];

% historyfiles         = cfg_files;
% historyfiles.name    = 'Monte Carlo history files';
% historyfiles.tag     = 'historyfiles';
% historyfiles.ufilter = {'.his','.mch'};    
% historyfiles.num     = [1 Inf];     
% historyfiles.help    = {'Select history files for this subject.'}; 

% Executable Branch
calculatePVE1      = cfg_exbranch;       
calculatePVE1.name = 'Calculate Partial Volume Effect';             
calculatePVE1.tag  = 'calculatePVE1';
calculatePVE1.val  = {NIRSmat DelPreviousData NewDirCopyNIRS dir_in}; 
calculatePVE1.prog = @nirs_run_calculatePVE;  
calculatePVE1.vout = @nirs_cfg_vout_calculatePVE; 
calculatePVE1.help = {'Calculate Partial Volume Effect'};

function vout = nirs_cfg_vout_calculatePVE(job)
vout = cfg_dep;                     
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reconstructions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

temp_pts       = cfg_entry; 
temp_pts.name    = 'Point temporel de l''inversion'; % The displayed name
temp_pts.tag     = 'temp_pts';       %file names
temp_pts.strtype = 'r';  
temp_pts.num     = [1 Inf];     % Number of inputs required 
% temp_pts.def = @(val)nirs_get_defaults('configMC1.nphotons', val{:});
temp_pts.help    = {'Input time point.'}; 
             
sens_vxsize= cfg_entry; 
sens_vxsize.name    = 'Voxel size in sensitivity matrix'; % The displayed name
sens_vxsize.tag     = 'sens_vxsize';       %file names
sens_vxsize.strtype = 'r';  
sens_vxsize.num     = [1 1];     % Number of inputs required 
% sens_vxsize.def = @(val)nirs_get_defaults('configMC1.nphotons', val{:});
sens_vxsize.help    = {'Input time point.'}; 

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

WLruns           = cfg_menu;
WLruns.name      = 'Runs';
WLruns.tag       = 'WLruns';
WLruns.labels    = {'One' 'Each WL separately'};
WLruns.values    = {1,2};
WLruns.val       = {1};
WLruns.help      = {''};

beta_wtd           = cfg_menu;
beta_wtd.name      = 'Beta : mua or Hbs';
beta_wtd.tag       = 'beta_wtd';
beta_wtd.labels    = {'mua' 'hbs'};
beta_wtd.values    = {1,2};
beta_wtd.val       = {1};
beta_wtd.help      = {'Choose Delta mua or Delta[HbO] and Delta[HBR]'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration: 3D reconstruction -- Tikhonov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha         = cfg_entry; 
alpha.name    = 'Hyperparameter';
alpha.tag     = 'alpha';
alpha.strtype = 'r';  
alpha.num     = [1 1];
alpha.val     = {1};
alpha.help    = {'Tunes the model :'
    'If small, the solution favors minimizing the residual error with the measured data.'
    'If large, the solution is biased towards matching the prior (no activity).'}';

tikh_method           = cfg_menu;
tikh_method.name      = 'Tikhonov regularization method';
tikh_method.tag       = 'tikh_method';
tikh_method.labels    = {'a l''ancienne !' 'Tikhonov' 'Extended Tikhonov' 'Simple Bayesian Interpretation'};
tikh_method.values    = {0,1,2,3};
tikh_method.val       = {1};
tikh_method.help      = {'Choose Tikhonov regularization reconstruction method (all taken from Hierarchical Bayesian regularization of reconstructions for diffuse optical tomography using multiple priors, Huppert).'
    '-- a l''ancienne is the first method working (false but kept to be able to compare results)'
    '-- Tikhonov is the simplest method of regularization'
    '-- Extended Tikhonov uses Li et al. model'
    '-- Simple Bayesian Interpretation uses covariances for the norms.'
    }';

% Executable Branch
tikhonov1      = cfg_exbranch;       
tikhonov1.name = 'Tikhonov inversion';             
tikhonov1.tag  = 'tikhonov1';
tikhonov1.val  = {NIRSmat NewDirCopyNIRS temp_pts dir_in sens_vxsize tikh_method alpha}; 
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
ReML_method           = cfg_menu;
ReML_method.name      = 'ReML method';
ReML_method.tag       = 'ReML_method';
ReML_method.labels    = {'Huppert' 'SPM'};
ReML_method.values    = {0,1};
ReML_method.val       = {0};
ReML_method.help      = {'Choose ReML reconstruction method.'};

%%%%%%%%%%%%%BOLD

% Executable Branch
ReMLreconstruct1      = cfg_exbranch;       
ReMLreconstruct1.name = '3D NIRS data ReML reconstruction';             
ReMLreconstruct1.tag  = 'ReMLreconstruct1';
ReMLreconstruct1.val  = {NIRSmat NewDirCopyNIRS beta_wtd temp_pts dir_in sens_vxsize ReML_method WLruns};   
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
% check reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outreconstruct_Hb         = cfg_files;
outreconstruct_Hb.name    = 'Select Hb output files';
outreconstruct_Hb.tag     = 'outreconstruct_Hb';
outreconstruct_Hb.ufilter = {'.nii'};    
outreconstruct_Hb.num     = [1 Inf];     
outreconstruct_Hb.help    = {'.'};

% Executable Branch
checkreconstruct1      = cfg_exbranch;       
checkreconstruct1.name = 'Check reconstruction';             
checkreconstruct1.tag  = 'checkreconstruct1';
checkreconstruct1.val  = {NIRSmat outreconstruct_Hb};   
checkreconstruct1.prog = @nirs_run_checkreconstruct;  
checkreconstruct1.vout = @nirs_cfg_vout_checkreconstruct; 
checkreconstruct1.help = {'Check reconstruction.'};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_checkreconstruct(job)
vout = cfg_dep;                     
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

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
% hrf Canonical HRF
% ---------------------------------------------------------------------
hrf         = cfg_branch;
hrf.tag     = 'hrf';
hrf.name    = 'Canonical HRF';
hrf.val     = {derivs };
hrf.help    = {'Canonical Hemodynamic Response Function. This is the default option. Contrasts of these effects have a physical interpretation and represent a parsimonious way of characterising event-related responses. This option is also useful if you wish to look separately at activations and deactivations (this is implemented using a t-contrast with a +1 or -1 entry over the canonical regressor). '};

% ---------------------------------------------------------------------
% length Window length
% ---------------------------------------------------------------------
length         = cfg_entry;
length.tag     = 'length';
length.name    = 'Window length';
length.help    = {'Post-stimulus window length (in seconds)'};
length.strtype = 'e';
length.num     = [1 1];
% ---------------------------------------------------------------------
% order Order
% ---------------------------------------------------------------------
order         = cfg_entry;
order.tag     = 'order';
order.name    = 'Order';
order.help    = {'Number of basis functions'};
order.strtype = 'e';
order.num     = [1 1];
% ---------------------------------------------------------------------
% gamma Gamma Functions
% ---------------------------------------------------------------------
gamma         = cfg_branch;
gamma.tag     = 'gamma';
gamma.name    = 'Gamma Functions';
gamma.val     = {length order };
gamma.help    = {'Gamma basis functions - requires SPM{F} for inference if more than one basis function.'};

% ---------------------------------------------------------------------
% bases Basis Functions
% ---------------------------------------------------------------------
bases         = cfg_choice;
bases.tag     = 'bases';
bases.name    = 'Basis Functions';
bases.val     = {hrf };
bases.help    = {'This option is only used for tests using gamma function. --- The most common choice of basis function is the Canonical HRF with or without time and dispersion derivatives. '};
bases.values  = {hrf gamma}; %{hrf fourier fourier_han gamma fir };

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
volt.val = {1};

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

hpf_butter_freq         = cfg_entry; 
hpf_butter_freq.name    = 'Cutoff frequency for HPF';
hpf_butter_freq.tag     = 'hpf_butter_freq';       
hpf_butter_freq.strtype = 'r';
hpf_butter_freq.num     = [1 1];     
hpf_butter_freq.def     = @(val)nirs_get_defaults(...
    'model_specify.wls_bglm_specify.hpf_butter.hpf_butter_On.hpf_butter_freq', val{:}); 
hpf_butter_freq.help    = {'Enter cutoff frequency in Hz for Butterworth HPF.'};

hpf_butter_order         = cfg_entry; 
hpf_butter_order.name    = 'Order of Butterworth HPF';
hpf_butter_order.tag     = 'hpf_butter_order';       
hpf_butter_order.strtype = 'r';
hpf_butter_order.num     = [1 1];     
hpf_butter_order.val     = {5};
hpf_butter_order.help    = {'Enter order of Butterworth HPF (preferred value = 3).'};

hpf_butter_On         = cfg_branch;
hpf_butter_On.tag     = 'hpf_butter_On';
hpf_butter_On.name    = 'Butterworth HP filter';
hpf_butter_On.val     = {hpf_butter_freq hpf_butter_order}; 
hpf_butter_On.help    = {'Butterworth high-pass filter.'};

hpf_butter_Off         = cfg_branch;
hpf_butter_Off.tag     = 'hpf_butter_Off';
hpf_butter_Off.name    = 'HP filter off';
hpf_butter_Off.val     = {}; 
hpf_butter_Off.help    = {'High pass filter turned off.'};

NoNIRSconfounds         = cfg_branch;
NoNIRSconfounds.tag     = 'NoNIRSconfounds';
NoNIRSconfounds.name    = 'No NIRS channels as confounds';
NoNIRSconfounds.val     = {}; 
NoNIRSconfounds.help    = {'No NIRS channels as confounds.'};

NumChConfounds         = cfg_entry; 
NumChConfounds.name    = 'Maximum Number of Confounds';
NumChConfounds.tag     = 'NumChConfounds';       
NumChConfounds.strtype = 'r';
NumChConfounds.num     = [1 1];     
NumChConfounds.val     = {1};
NumChConfounds.help    = {'Enter maximum number of NIRS channels to be included as physiological confounds.'};

MinChDist         = cfg_entry; 
MinChDist.name    = 'Minimum Channel Length';
MinChDist.tag     = 'MinChDist';       
MinChDist.strtype = 'r';
MinChDist.num     = [1 1];     
MinChDist.val     = {1.5};
MinChDist.help    = {'Enter minimum channel length allowed for inclusion as confound.'};

MaxChDist         = cfg_entry; 
MaxChDist.name    = 'maximum Channel Length';
MaxChDist.tag     = 'MaxChDist';       
MaxChDist.strtype = 'r';
MaxChDist.num     = [1 1];     
MaxChDist.val     = {2.5};
MaxChDist.help    = {'Enter maximum channel length allowed for inclusion as confound.'};

NumChConfounds         = cfg_entry; 
NumChConfounds.name    = 'Maximum Number of Confounds';
NumChConfounds.tag     = 'NumChConfounds';       
NumChConfounds.strtype = 'r';
NumChConfounds.num     = [1 1];     
NumChConfounds.val     = {1};
NumChConfounds.help    = {'Enter maximum number of NIRS channels to be included as physiological confounds.'};

NIRSconfounds         = cfg_branch;
NIRSconfounds.tag     = 'NIRSconfounds';
NIRSconfounds.name    = 'NIRS channels as confounds';
NIRSconfounds.val     = {NumChConfounds MinChDist MaxChDist}; 
NIRSconfounds.help    = {'NIRS channels as confounds.'};

NIRSchannelsConfound         = cfg_choice;
NIRSchannelsConfound.tag     = 'NIRSchannelsConfound';
NIRSchannelsConfound.name    = 'NIRS channels as confounds';
NIRSchannelsConfound.values  = {NoNIRSconfounds NIRSconfounds}; 
NIRSchannelsConfound.val     = {NoNIRSconfounds}; 
NIRSchannelsConfound.help    = {'NIRS channels  as confound regressors.'
    'When using this option, selected channels will be filtered with the'
    'Same parameters as for the GLM and included as confound regressors'
    'This may be useful to remove physiological noise'
    'Only HbO channels will be used'}';

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
filter_design_matrix.val = {1};
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
wls_bglm_specify.val  = {NIRSmat dir1 subj units time_res derivs bases ...
    volt GLM_include_cardiac GLM_include_Mayer NIRSchannelsConfound GenerateHbT flag_window ...
    channel_pca hpf_butter generate_trRV filter_design_matrix ...
     wls_or_bglm LiomDeleteLarge}; 
wls_bglm_specify.prog = @nirs_run_liom_GLM_specify;  
wls_bglm_specify.vout = @nirs_cfg_vout_liom_GLM_specify; 
wls_bglm_specify.help = {'Specify LIOM General Linear Model.'};

function vout = nirs_cfg_vout_liom_GLM_specify(job)
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
wls_bglm_estimate.prog = @nirs_run_liom_GLM_estimate;  
wls_bglm_estimate.vout = @nirs_cfg_vout_liom_GLM_estimate; 
wls_bglm_estimate.help = {'LIOM GLM Estimation: WLS (wavelet least square)'
            'and Bayesian GLM.'}';

function vout = nirs_cfg_vout_liom_GLM_estimate(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%Select view
view         = cfg_entry; 
view.name    = 'View'; 
view.tag     = 'view';    
view.strtype = 'r'; 
view.num     = [1 Inf];  
view.val     = {5};
view.help    = {['Enter view.  ',...
    '1: ventral  ',...
    '2: dorsal  ',...
    '3: right  ',...
    '4: left  ',...
    '5: frontal  ',...
    '6: occipital']}; % help text displayed


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % NIRS_SPM Contrast calculations
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%File for channel coregistration 'preproc_info' generated by NIRS_SPM
NIRS_SPM_Coregistration_Channels         = cfg_files;  
NIRS_SPM_Coregistration_Channels.name    = 'Select file of coregistration info'; 
NIRS_SPM_Coregistration_Channels.tag     = 'NIRS_SPM_Coregistration_Channels';       
NIRS_SPM_Coregistration_Channels.filter  = 'mat';    
NIRS_SPM_Coregistration_Channels.num     = [1 1];    
NIRS_SPM_Coregistration_Channels.help    = {'Select file of channel '
    'coregistration ''preproc_info'' generated by NIRS_SPM.'}'; 

% 
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
% 
% contrast_struct         = cfg_repeat;
% contrast_struct.tag     = 'contrast_struct';
% contrast_struct.name    = 'Contrasts';
% contrast_struct.help    = {'Specify contrasts'};
% contrast_struct.values  = {contrast_data};
% contrast_struct.num     = [1 Inf];
% 
% % Executable Branch
% NIRS_SPM_contrast      = cfg_exbranch;      
% NIRS_SPM_contrast.name = 'NIRS_SPM Contrast Calculations';            
% NIRS_SPM_contrast.tag  = 'NIRS_SPM_contrast';
% NIRS_SPM_contrast.val  = {Dmx_files NIRS_SPM_Coregistration_Channels ...
%                 view contrast_struct}; 
% NIRS_SPM_contrast.prog = @nirs_run_NIRS_SPM_contrast;  
% NIRS_SPM_contrast.vout = @nirs_cfg_vout_NIRS_SPM_contrast; 
% NIRS_SPM_contrast.help = {'NIRS_SPM Contrast Calculations.'};
% 
% function vout = nirs_cfg_vout_NIRS_SPM_contrast(job)
% vout = cfg_dep;                     
% vout.sname      = 'NIRS.mat';       
% vout.src_output = substruct('.','NIRSmat'); 
% vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Liom Contrast calculations - based on tube formula and code by NIRS_SPM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% liom_contrast_struct         = cfg_repeat;
% liom_contrast_struct.tag     = 'liom_contrast_struct';
% liom_contrast_struct.name    = 'Contrasts';
% liom_contrast_struct.help    = {'Specify contrasts'
%     'NOTE: this is an older method to specify contrasts, which will disappear.'}';
% liom_contrast_struct.values  = {contrast_data};
% liom_contrast_struct.num     = [0 Inf];

contrast_p_value         = cfg_entry; 
contrast_p_value.name    = 'Contrast uncorrected p_value';
contrast_p_value.tag     = 'contrast_p_value';       
contrast_p_value.strtype = 'r';
contrast_p_value.num     = [1 1];
contrast_p_value.val     = {0.05};
contrast_p_value.help    = {'Contrast uncorrected p_value'}; 

spatial_LPF_radius         = cfg_entry; 
spatial_LPF_radius.name    = 'Spatial LPF radius';
spatial_LPF_radius.tag     = 'spatial_LPF_radius';       
spatial_LPF_radius.strtype = 'r';
spatial_LPF_radius.num     = [1 1];     
spatial_LPF_radius.val     = {3}; 
spatial_LPF_radius.help    = {'Enter radius of spatial low pass filter in pixels.'
    'One pixel is very approximately 1 mm. FWHM will be twice this radius.'
    'If 0 is entered, this is equivalent to no spatial filtering.'
    'Spatial filtering will be applied linearly even though '
    'stereographic projection is nonlinear.'}';

spatial_LPF_On         = cfg_branch;
spatial_LPF_On.tag     = 'spatial_LPF_On';
spatial_LPF_On.name    = 'Spatial LP filter';
spatial_LPF_On.val     = {spatial_LPF_radius}; 
spatial_LPF_On.help    = {'Spatial low-pass filter.'};

spatial_LPF_Off         = cfg_branch;
spatial_LPF_Off.tag     = 'spatial_LPF_Off';
spatial_LPF_Off.name    = 'Spatial filter off';
spatial_LPF_Off.val     = {}; 
spatial_LPF_Off.help    = {'Spatial low pass filter turned off.'};

spatial_LPF      = cfg_choice;
spatial_LPF.tag  = 'spatial_LPF';
spatial_LPF.name = 'Spatial Low Pass Filter';
spatial_LPF.values = {spatial_LPF_On spatial_LPF_Off};
spatial_LPF.val = {spatial_LPF_Off};
spatial_LPF.help = {'Choose whether to include a spatial Low Pass Filter'
    'on the interpolated estimates and their variance before constructing'
    'statistical maps.'}';

contrast_figures      = cfg_menu;
contrast_figures.tag  = 'contrast_figures';
contrast_figures.name = 'Generate figures';
contrast_figures.labels = {'No','Both .fig and .tiff','Only .fig','Only .tiff'};
contrast_figures.values = {0,1,2,3};
contrast_figures.val = {3};
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


colorbar_max2         = cfg_entry; 
colorbar_max2.name    = 'Colorbar maximum value';
colorbar_max2.tag     = 'colorbar_max2';       
colorbar_max2.strtype = 'r';
colorbar_max2.num     = [1 1];
colorbar_max2.val     = {-2};
colorbar_max2.help    = {'Enter maximum value for colorbar for negative maps'}; 

colorbar_min2         = cfg_entry; 
colorbar_min2.name    = 'Colorbar minimum value';
colorbar_min2.tag     = 'colorbar_min2';       
colorbar_min2.strtype = 'r';
colorbar_min2.num     = [1 1];
colorbar_min2.val     = {-5};
colorbar_min2.help    = {'Enter minimum value for colorbar for negative maps'}; 

colorbar_override      = cfg_branch;
colorbar_override.name      = 'Override colorbar';
colorbar_override.tag       = 'colorbar_override';
colorbar_override.val       = {colorbar_min colorbar_max colorbar_min2 colorbar_max2}; 
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

GroupColorbars      = cfg_menu;
GroupColorbars.tag  = 'GroupColorbars';
GroupColorbars.name = 'Group or separate colorbars';
GroupColorbars.labels = {'Group','Separate'};
GroupColorbars.values = {1,0};
GroupColorbars.val = {0};
GroupColorbars.help = {'When considering deactivations (inverted responses),'
    'This allows displaying one (group) or two (separate) colorbars'}';

SmallFigures      = cfg_menu;
SmallFigures.tag  = 'SmallFigures';
SmallFigures.name = 'Large or small figures';
SmallFigures.labels = {'Large','Small'};
SmallFigures.values = {0,1};
SmallFigures.val = {1};
SmallFigures.help = {'Write to disk large or small (compressed) figures.'}';

ProcessContrastsBySession      = cfg_menu;
ProcessContrastsBySession.tag  = 'ProcessContrastsBySession';
ProcessContrastsBySession.name = 'Process Contrasts By Session';
ProcessContrastsBySession.labels = {'Yes','No','Both'};
ProcessContrastsBySession.values = {1,0,2};
ProcessContrastsBySession.val = {1};
ProcessContrastsBySession.help = {'Important option, careful; only applies'
    'when there are two or more sessions'
    'if yes, contrasts defined over more than one session will be ignored'
    'if no, contrasts defined over only one session will be processed as'
    'a contrast defined over the full design matrix over all sessions, with'
    'statistics, e.g. number of degrees of freedom, calculated for the full design matrix.'
    'if both, both options will be run.'}';

GroupMultiSession           = cfg_menu;
GroupMultiSession.name      = 'Group Multi-Session';
GroupMultiSession.tag       = 'GroupMultiSession';
GroupMultiSession.labels    = {'Yes' 'No'};
GroupMultiSession.values    = {1,0};
GroupMultiSession.val       = {0};
GroupMultiSession.help      = {'Group Multi Session'
    'If selected, with option ProcessContrastBySession set to 0,'
    'Contrasts defined as vectors over all sessions will be treated.'}';

% Study_type           = cfg_menu;
% Study_type.name      = 'Study type';
% Study_type.tag       = 'Study_type';
% Study_type.labels    = {'Group Single Session' 'Group Multi Sessions' 'Single Subject Multi Sessions'};
% Study_type.values    = {2,1,0};
% Study_type.val       = {2};
% Study_type.help      = {'Study type: note that all cases can be used for one or more subjects.'
%     'However, the treatment by liom_contrast and/or liom_group will be affected by this choice'
%     'Group Single Session: use when there is only one session, or with more than one session but contrasts do not combine different sessions'
%     'Group Multi Sessions: use when contrasts are combined over more than one session'
%     'Single Subject Multi Sessions: use when subjects will not be combined into a group analysis.'}';

GroupFiguresIntoSubplots      = cfg_menu;
GroupFiguresIntoSubplots.tag  = 'GroupFiguresIntoSubplots';
GroupFiguresIntoSubplots.name = 'Group Figures Into Subplots';
GroupFiguresIntoSubplots.labels = {'Yes','No'};
GroupFiguresIntoSubplots.values = {1,0};
GroupFiguresIntoSubplots.val = {1};
GroupFiguresIntoSubplots.help = {'Group Figures Into Subplots.'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From spm spm_cfg_con
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------------------------------------------------------------------
% name Name
% ---------------------------------------------------------------------
name         = cfg_entry;
name.tag     = 'name';
name.name    = 'Name';
name.help    = {'Name of contrast'};
name.strtype = 's';
name.num     = [1 Inf];
% ---------------------------------------------------------------------
% convec T contrast vector
% ---------------------------------------------------------------------
convec         = cfg_entry;
convec.tag     = 'convec';
convec.name    = 'T contrast vector';
convec.help    = {'Enter T contrast vector. This is done similarly to the contrast manager. A 1 x n vector should be entered for T-contrasts.'};
convec.strtype = 'e';
convec.num     = [1 Inf];
% ---------------------------------------------------------------------
% sessrep Replicate over sessions
% ---------------------------------------------------------------------
sessrep         = cfg_menu;
sessrep.tag     = 'sessrep';
sessrep.name    = 'Replicate over sessions';
sessrep.val = {'none'};
sessrep.help    = {
                   'If there are multiple sessions with identical conditions, one might want to specify contrasts which are identical over sessions. This can be done automatically based on the contrast spec for one session.'
                   'Contrasts can be either replicated (thus testing average effects over sessions) or created per session. In both cases, zero padding up to the length of each session and the block effects is done automatically. In addition, weights of replicated contrasts can be scaled by the number of sessions. This allows to use the same contrast manager batch for fMRI analyses with a variable number of sessions. The scaled contrasts can then be compared in a 2nd level model without a need for further adjustment of effect sizes.'
}';
sessrep.labels = {
                  'Don''t replicate'
                  'Replicate'
                  'Replicate&Scale'
                  'Create per session'
                  'Both: Replicate + Create per session'
                  'Both: Replicate&Scale + Create per session'
}';
sessrep.values = {
                  'none'
                  'repl'
                  'replsc'
                  'sess'
                  'both'
                  'bothsc'
}';
% ---------------------------------------------------------------------
% tcon T-contrast
% ---------------------------------------------------------------------
tcon         = cfg_branch;
tcon.tag     = 'tcon';
tcon.name    = 'T-contrast';
tcon.val     = {name convec sessrep };
tcon.help    = {
                '* Simple one-dimensional contrasts for an SPM{T}'
                ''
                'A simple contrast for an SPM{T} tests the null hypothesis c''B=0 against the one-sided alternative c''B>0, where c is a column vector. '
                ''
                '    Note that throughout SPM, the transpose of the contrast weights is used for display and input. That is, you''ll enter and visualise c''. For an SPM{T} this will be a row vector.'
                ''
                'For example, if you have a design in which the first two columns of the design matrix correspond to the effects for "baseline" and "active" conditions respectively, then a contrast with weights c''=[-1,+1,0,...] (with zero weights for any other parameters) tests the hypothesis that there is no "activation" (the parameters for both conditions are the same), against the alternative that there is some activation (i.e. the parameter for the "active" condition is greater than that for the "baseline" condition). The resulting SPM{T} (created by spm_getSPM.m) is a statistic image, with voxel values the value of the t-statistic for the specified contrast at that location. Areas of the SPM{T} with high voxel values indicate evidence for "activation". To look for areas of relative "de-activation", the inverse contrast could be used c''=[+1,-1,0,...].'
                ''
                'Similarly, if you have a design where the third column in the design matrix is a covariate, then the corresponding parameter is essentially a regression slope, and a contrast with weights c''=[0,0,1,0,...] (with zero weights for all parameters but the third) tests the hypothesis of zero regression slope, against the alternative of a positive slope. This is equivalent to a test no correlation, against the alternative of positive correlation. If there are other terms in the model beyond a constant term and the covariate, then this correlation is apartial correlation, the correlation between the data Y and the covariate, after accounting for the other effects.'
}';
% ---------------------------------------------------------------------
% name Name
% ---------------------------------------------------------------------
name         = cfg_entry;
name.tag     = 'name';
name.name    = 'Name';
name.help    = {'Name of contrast'};
name.strtype = 's';
name.num     = [1 Inf];
% ---------------------------------------------------------------------
% convec F contrast vector
% ---------------------------------------------------------------------
convec         = cfg_entry;
convec.tag     = 'convec';
convec.name    = 'F contrast vector';
convec.help    = {'Enter F contrast vector. This is done similarly to the contrast manager. One or multiline contrasts may be entered.'};
convec.strtype = 'e';
convec.num     = [Inf Inf];
% ---------------------------------------------------------------------
% generic Contrast vectors
% ---------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Contrast vectors';
generic.help    = {'F contrasts are defined by a series of vectors.'};
generic.values  = {convec };
generic.num     = [1 Inf];
% ---------------------------------------------------------------------
% sessrep Replicate over sessions
% ---------------------------------------------------------------------
sessrep         = cfg_menu;
sessrep.tag     = 'sessrep';
sessrep.name    = 'Replicate over sessions';
sessrep.val = {'none'};
sessrep.help    = {
                   'If there are multiple sessions with identical conditions, one might want to specify contrasts which are identical over sessions. This can be done automatically based on the contrast spec for one session.'
                   'Contrasts can be either replicated (either testing average effects over sessions or per-session/condition effects) or created per session. In both cases, zero padding up to the length of each session and the block effects is done automatically.'
}';
sessrep.labels = {
                  'Don''t replicate'
                  'Replicate (average over sessions)'
                  'Replicate (no averaging)'
                  'Create per session'
                  'Both - ''Per session'' and ''Replicate (average over sessions)'''
}';
sessrep.values = {
                  'none'
                  'repl'
                  'replna'
                  'sess'
                  'both'
}';
% ---------------------------------------------------------------------
% fcon F-contrast
% ---------------------------------------------------------------------
fcon         = cfg_branch;
fcon.tag     = 'fcon';
fcon.name    = 'F-contrast';
fcon.val     = {name generic sessrep };
fcon.help    = {
                '* Linear constraining matrices for an SPM{F}'
                ''
                'The null hypothesis c''B=0 can be thought of as a (linear) constraint on the full model under consideration, yielding a reduced model. Taken from the viewpoint of two designs, with the full model an extension of the reduced model, the null hypothesis is that the additional terms in the full model are redundent.'
                ''
                'Statistical inference proceeds by comparing the additional variance explained by full design over and above the reduced design to the error variance (of the full design), an "Extra Sum-of-Squares" approach yielding an F-statistic for each voxel, whence an SPM{F}.'
                ''
                'This is useful in a number of situations:'
                ''
                '* Two sided tests'
                ''
                'The simplest use of F-contrasts is to effect a two-sided test of a simple linear contrast c''B, where c is a column vector. The SPM{F} is the square of the corresponding SPM{T}. High values of the SPM{F} therefore indicate evidence against the null hypothesis c''B=0 in favour of the two-sided alternative c''B~=0.'
                ''
                '* General linear hypotheses'
                ''
                'Where the contrast weights is a matrix, the rows of the (transposed) contrast weights matrix c'' must define contrasts in their own right, and the test is effectively simultaneously testing the null hypotheses associated with the individual component contrasts with weights defined in the rows. The null hypothesis is still c''B=0, but since c is a matrix, 0 here is a zero vector rather than a scalar zero, asserting that under the null hypothesis all the component hypotheses are true.'
                ''
                'For example: Suppose you have a language study with 3 word categories (A,B & C), and would like to test whether there is any difference at all between the three levels of the "word category" factor.'
                ''
                'The design matrix might look something like:'
                ''
                '         [ 1 0 0 ..]'
                '         [ : : : ..]'
                '         [ 1 0 0 ..]'
                '         [ 0 1 0 ..]'
                '    X =  [ : : : ..]'
                '         [ 0 1 0 ..]'
                '         [ 0 0 1 ..]'
                '         [ : : : ..]'
                '         [ 0 0 1 ..]'
                '         [ 0 0 0 ..]'
                '         [ : : : ..]'
                ''
                ' ...with the three levels of the "word category" factor modelled in the  first three columns of the design matrix.'
                ''
                'The matrix of contrast weights will look like:'
                ''
                ' c'' = [1 -1  0 ...;'
                '       0  1 -1 ...]'
                ''
                'Reading the contrasts weights in each row of c'', we see that row 1 states that category A elicits the same response as category B, row 2 that category B elicits the same response as category C, and hence together than categories A, B & C all elicit the same response.'
                ''
                'The alternative hypothesis is simply that the three levels are not all the same, i.e. that there is some difference in the paraeters for the three levels of the factor: The first and the second categories produce different brain responses, OR the second and third categories, or both.'
                ''
                'In other words, under the null hypothesis (the categories produce the same brain responses), the model reduces to one in which the three level "word category" factor can be replaced by a single "word" effect, since there is no difference in the parameters for each category. The corresponding design matrix would have the first three columns replaced by a single column that is the sum (across rows) of the first three columns in the design matric above, modelling the brain response to a word, whatever is the category. The F-contrast above is in fact testing the hypothesis that this reduced design doesn''t account for significantly less variance than the full design with an effect for each word category.'
                ''
                'Another way of seeing that, is to consider a reparameterisation of the model, where the first column models effects common to all three categories, with the second and third columns modelling the differences between the three conditions, for example:'
                ''
                '         [ 1  1  0 ..]'
                '         [ :  :  : ..]'
                '         [ 1  1  0 ..]'
                '         [ 1  0  1 ..]'
                '    X =  [ :  :  : ..]'
                '         [ 1  0  1 ..]'
                '         [ 1 -1 -1 ..]'
                '         [ :  :  : ..]'
                '         [ 1 -1 -1 ..]'
                '         [ 0  0  0 ..]'
                '         [ :  :  : ..]'
                ''
                'In this case, an equivalent F contrast is of the form'
                ' c'' = [ 0 1 0 ...;'
                '        0 0 1 ...]'
                'and would be exactly equivalent to the previous contrast applied to the previous design. In this latter formulation, you are asking whewher the two columns modelling the "interaction space" account for a significant amount of variation (variance) of the data. Here the component contrasts in the rows of c'' are simply specifying that the parameters for the corresponding rows are are zero, and it is clear that the F-test is comparing this full model with a reduced model in which the second and third columns of X are omitted.'
                ''
                '    Note the difference between the following two F-contrasts:'
                '         c'' = [ 0 1 0 ...;     (1)'
                '                0 0 1 ...]'
                '     and'
                '         c'' = [ 0 1 1 ...]     (2)'
                ''
                '    The first is an F-contrast, testing whether either of the parameters for the effects modelled in the 2nd & 3rd columns of the design matrix are significantly different from zero. Under the null hypothesis c''B=0, the first contrast imposes a two-dimensional constraint on the design. The second contrast tests whether the SUM of the parameters for the 2nd & 3rd columns is significantly different from zero. Under the null hypothesis c''B=0, this second contrast only imposes a one dimensional constraint on the design.'
                ''
                '    An example of the difference between the two is that the first contrast would be sensitive to the situation where the 2nd & 3rd parameters were +a and -a, for some constant a, wheras the second contrast would not detect this, since the parameters sum to zero.'
                ''
                'The test for an effect of the factor "word category" is an F-test with 3-1=2 "dimensions", or degrees of freedom.'
                ''
                '* Testing the significance of effects modelled by multiple columns'
                ''
                'A conceptially similar situation arises when one wonders whether a set of coufound effects are explaining any variance in the data. One important advantage of testing the with F contrasts rather than one by one using SPM{T}''s is the following. Say you have two covariates that you would like to know whether they can "predict" the brain responses, and these two are correlated (even a small correlation would be important in this instance). Testing one and then the other may lead you to conclude that there is no effect. However, testing with an F test the two covariates may very well show a not suspected effect. This is because by testing one covariate after the other, one never tests for what is COMMON to these covariates (see Andrade et al, Ambiguous results in functional neuroimaging, NeuroImage, 1999).'
                ''
                ''
                'More generally, F-tests reflect the usual analysis of variance, while t-tests are traditionally post hoc tests, useful to see in which direction is an effect going (positive or negative). The introduction of F-tests can also be viewed as a first means to do model selection.'
                ''
                ''
                'Technically speaking, an F-contrast defines a number of directions (as many as the rank of the contrast) in the space spanned by the column vectors of the design matrix. These directions are simply given by X*c if the vectors of X are orthogonal, if not, the space define by c is a bit more complex and takes care of the correlation within the design matrix. In essence, an F-contrast is defining a reduced model by imposing some linear constraints (that have to be estimable, see below) on the parameters estimates. Sometimes, this reduced model is simply made of a subset of the column of the original design matrix but generally, it is defined by a combination of those columns. (see spm_FcUtil for what (I hope) is an efficient handling of F-contrats computation).'
}';
% ---------------------------------------------------------------------
% name Name
% ---------------------------------------------------------------------
name         = cfg_entry;
name.tag     = 'name';
name.name    = 'Name';
name.help    = {'Name of contrast'};
name.strtype = 's';
name.num     = [1 Inf];
% ---------------------------------------------------------------------
% conweight Contrast weight
% ---------------------------------------------------------------------
conweight         = cfg_entry;
conweight.tag     = 'conweight';
conweight.name    = 'Contrast weight';
conweight.help    = {'The contrast weight for the selected column.'};
conweight.strtype = 'e';
conweight.num     = [1 1];
% ---------------------------------------------------------------------
% colcond Condition #
% ---------------------------------------------------------------------
colcond         = cfg_entry;
colcond.tag     = 'colcond';
colcond.name    = 'Condition #';
colcond.help    = {'Select which condition function set is to be contrasted.'};
colcond.strtype = 'e';
colcond.num     = [1 1];
% ---------------------------------------------------------------------
% colbf Basis function #
% ---------------------------------------------------------------------
colbf         = cfg_entry;
colbf.tag     = 'colbf';
colbf.name    = 'Basis function #';
colbf.help    = {'Select which basis function from the basis function set is to be contrasted.'};
colbf.strtype = 'e';
colbf.num     = [1 1];
% ---------------------------------------------------------------------
% colmod Parametric modulation #
% ---------------------------------------------------------------------
colmod         = cfg_entry;
colmod.tag     = 'colmod';
colmod.name    = 'Parametric modulation #';
colmod.help    = {'Select which parametric modulation is to be contrasted. If there is no time/parametric modulation, enter "1". If there are both time and parametric modulations, then time modulation comes before parametric modulation.'};
colmod.strtype = 'e';
colmod.num     = [1 1];
% ---------------------------------------------------------------------
% colmodord Parametric modulation order
% ---------------------------------------------------------------------
colmodord         = cfg_entry;
colmodord.tag     = 'colmodord';
colmodord.name    = 'Parametric modulation order';
colmodord.help    = {
                     'Order of parametric modulation to be contrasted. '
                     ''
                     '0 - the basis function itself, 1 - 1st order mod etc'
}';
colmodord.strtype = 'e';
colmodord.num     = [1 1];
% ---------------------------------------------------------------------
% colconds Contrast entry
% ---------------------------------------------------------------------
colconds         = cfg_branch;
colconds.tag     = 'colconds';
colconds.name    = 'Contrast entry';
colconds.val     = {conweight colcond colbf colmod colmodord };
colconds.help    = {''};
% ---------------------------------------------------------------------
% generic T contrast for conditions
% ---------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'T contrast for conditions';
generic.help    = {'Assemble your contrast column by column.'};
generic.values  = {colconds };
generic.num     = [1 Inf];
% ---------------------------------------------------------------------
% colreg T contrast for extra regressors
% ---------------------------------------------------------------------
colreg         = cfg_entry;
colreg.tag     = 'colreg';
colreg.name    = 'T contrast for extra regressors';
colreg.help    = {'Enter T contrast vector for extra regressors.'};
colreg.strtype = 'e';
colreg.num     = [1 Inf];
% ---------------------------------------------------------------------
% coltype Contrast columns
% ---------------------------------------------------------------------
coltype         = cfg_choice;
coltype.tag     = 'coltype';
coltype.name    = 'Contrast columns';
coltype.val     = {generic };
coltype.help    = {'Contrasts can be specified either over conditions or over extra regressors.'};
coltype.values  = {generic colreg };
% ---------------------------------------------------------------------
% sessions Session(s)
% ---------------------------------------------------------------------
sessions         = cfg_entry;
sessions.tag     = 'sessions';
sessions.name    = 'Session(s)';
sessions.help    = {'Enter session number(s) for which this contrast should be created. If more than one session number is specified, the contrast will be an average contrast over the specified conditions or regressors from these sessions.'};
sessions.strtype = 'e';
sessions.num     = [1 Inf];
% ---------------------------------------------------------------------
% tconsess T-contrast (cond/sess based)
% ---------------------------------------------------------------------
tconsess         = cfg_branch;
tconsess.tag     = 'tconsess';
tconsess.name    = 'T-contrast (cond/sess based)';
tconsess.val     = {name coltype sessions };
tconsess.help    = {
                    'Define a contrast in terms of conditions or regressors instead of columns of the design matrix. This allows to create contrasts automatically even if some columns are not always present (e.g. parametric modulations).'
                    ''
                    'Each contrast column can be addressed by specifying'
                    '* session number'
                    '* condition number'
                    '* basis function number'
                    '* parametric modulation number and'
                    '* parametric modulation order.'
                    ''
                    'If the design is specified without time or parametric modulation, SPM creates a "pseudo-modulation" with order zero. To put a contrast weight on a basis function one therefore has to enter "1" for parametric modulation number and "0" for parametric modulation order.'
                    ''
                    'Time and parametric modulations are not distinguished internally. If time modulation is present, it will be parametric modulation "1", and additional parametric modulations will be numbered starting with "2".'
                    ''
                    '* Simple one-dimensional contrasts for an SPM{T}'
                    ''
                    'A simple contrast for an SPM{T} tests the null hypothesis c''B=0 against the one-sided alternative c''B>0, where c is a column vector. '
                    ''
                    '    Note that throughout SPM, the transpose of the contrast weights is used for display and input. That is, you''ll enter and visualise c''. For an SPM{T} this will be a row vector.'
                    ''
                    'For example, if you have a design in which the first two columns of the design matrix correspond to the effects for "baseline" and "active" conditions respectively, then a contrast with weights c''=[-1,+1,0,...] (with zero weights for any other parameters) tests the hypothesis that there is no "activation" (the parameters for both conditions are the same), against the alternative that there is some activation (i.e. the parameter for the "active" condition is greater than that for the "baseline" condition). The resulting SPM{T} (created by spm_getSPM.m) is a statistic image, with voxel values the value of the t-statistic for the specified contrast at that location. Areas of the SPM{T} with high voxel values indicate evidence for "activation". To look for areas of relative "de-activation", the inverse contrast could be used c''=[+1,-1,0,...].'
                    ''
                    'Similarly, if you have a design where the third column in the design matrix is a covariate, then the corresponding parameter is essentially a regression slope, and a contrast with weights c''=[0,0,1,0,...] (with zero weights for all parameters but the third) tests the hypothesis of zero regression slope, against the alternative of a positive slope. This is equivalent to a test no correlation, against the alternative of positive correlation. If there are other terms in the model beyond a constant term and the covariate, then this correlation is apartial correlation, the correlation between the data Y and the covariate, after accounting for the other effects.'
}';
% ---------------------------------------------------------------------
% consess Contrast Sessions
% ---------------------------------------------------------------------
consess         = cfg_repeat;
consess.tag     = 'consess';
consess.name    = 'Contrast Sessions';
consess.help    = { ' NOTE: IF ONLY BASIC CONTRASTS ARE REQUIRED,'
                    ' No need to specify any contrasts -- they will be generated automatically.'
                    ' If more contrasts are required, all need to be specified.'
                    ' ' 
                   'For general linear model Y = XB + E with data Y, desgin matrix X, parameter vector B, and (independent) errors E, a contrast is a linear combination of the parameters c''B. Usually c is a column vector, defining a simple contrast of the parameters, assessed via an SPM{T}. More generally, c can be a matrix (a linear constraining matrix), defining an "F-contrast" assessed via an SPM{F}.'
                   ''
                   'The vector/matrix c contains the contrast weights. It is this contrast weights vector/matrix that must be specified to define the contrast. The null hypothesis is that the linear combination c''B is zero. The order of the parameters in the parameter (column) vector B, and hence the order to which parameters are referenced in the contrast weights vector c, is determined by the construction of the design matrix.'
                   ''
                   'There are two types of contrast in SPM: simple contrasts for SPM{T}, and "F-contrasts" for SPM{F}.'
                   ''
                   'For a thorough theoretical treatment, see the Human Brain Function book and the statistical literature referenced therein.'
                   ''
                   ''
                   '* Non-orthogonal designs'
                   ''
                   'Note that parameters zero-weighted in the contrast are still included in the model. This is particularly important if the design is not orthogonal (i.e. the columns of the design matrix are not orthogonal). In effect, the significance of the contrast is assessed *after* accounting for the other effects in the design matrix. Thus, if two covariates are correlated, testing the significance of the parameter associated with one will only test for the part that is not present in the second covariate. This is a general point that is also true for F-contrasts. See Andrade et al, Ambiguous results in functional neuroimaging, NeuroImage, 1999, for a full description of the effect of non othogonal design testing.'
                   ''
                   ''
                   '* Estimability'
                   ''
                   'The contrast c''B is estimated by c''b, where b are the parameter estimates given by b=pinv(X)*Y.'
                   ''
                   'However, if a design is rank-deficient (i.e. the columns of the design matrix are not linearly independent), then the parameters are not unique, and not all linear combinations of the parameter are valid contrasts, since contrasts must be uniquely estimable.'
                   ''
                   'A weights vector defines a valid contrast if and only if it can be constructed as a linear combination of the rows of the design matrix. That is c'' (the transposed contrast vector - a row vector) is in the row-space of the design matrix.'
                   ''
                   'Usually, a valid contrast will have weights that sum to zero over the levels of a factor (such as condition).'
                   ''
                   'A simple example is a simple two condition design including a constant, with design matrix'
                   ''
                   '          [ 1 0 1 ]'
                   '          [ : : : ]'
                   '     X =  [ 1 0 1 ]'
                   '          [ 0 1 1 ]'
                   '          [ : : : ]'
                   '          [ 0 1 1 ]'
                   ''
                   'The first column corresponds to condition 1, the second to condition 2, and the third to a constant (mean) term. Although there are three columns to the design matrix, the design only has two degrees of freedom, since any one column can be derived from the other two (for instance, the third column is the sum of the first two). There is no unique set of parameters for this model, since for any set of parameters adding a constant to the two condition effects and subtracting it from the constant effect yields another set of viable parameters. However, the difference between the two condition effects is uniquely estimated, so c''=[-1,+1,0] does define a contrast.'
                   ''
                   'If a parameter is estimable, then the weights vector with a single "1" corresponding to that parameter (and zero elsewhere) defines a valid contrast.'
                   ''
                   ''
                   '* Multiple comparisons'
                   ''
                   'Note that SPM implements no corrections to account for you looking at multiple contrasts.'
                   ''
                   'If you are interested in a set of hypotheses that together define a consistent question, then you should account for this when assessing the individual contrasts. A simple Bonferroni approach would assess N simultaneous contrasts at significance level alpha/N, where alpha is the chosen significance level (usually 0.05).'
                   ''
                   'For two sided t-tests using SPM{T}s, the significance level should be halved. When considering both SPM{T}s produced by a contrast and it''s inverse (the contrast with negative weights), to effect a two-sided test to look for both "increases" and "decreases", you should review each SPM{T} at at level 0.05/2 rather than 0.05. (Or consider an F-contrast!)'
                   ''
                   ''
                   '* Contrast images and ESS images'
                   ''
                   'For a simple contrast, SPM (spm_getSPM.m) writes a contrast image: con_????.{img,nii}, with voxel values c''b. (The ???? in the image names are replaced with the contrast number.) These contrast images (for appropriate contrasts) are suitable summary images of an effect at this level, and can be used as input at a higher level when effecting a random effects analysis. See spm_RandFX.man for further details.'
                   ''
                   'For an F-contrast, SPM (spm_getSPM.m) writes the Extra Sum-of-Squares (the difference in the residual sums of squares for the full and reduced model) as ess_????.{img,nii}. (Note that the ess_????.{img,nii} and SPM{T,F}_????.{img,nii} images are not suitable input for a higher level analysis.)'
}';
consess.values  = {tcon fcon tconsess };
consess.num     = [0 Inf];

%%%%%%%%%%%%
% end - from spm_cfg_con 
%%%%%%%%%%%
output_unc      = cfg_menu;
output_unc.tag  = 'output_unc';
output_unc.name = 'Output unc. figures';
output_unc.labels = {'Yes','No'};
output_unc.values = {1,0};
output_unc.val = {0};
output_unc.help = {'Output figures that are uncorrected against false positives.'}';

write_neg_pos      = cfg_menu;
write_neg_pos.tag  = 'write_neg_pos';
write_neg_pos.name = 'Write neg/pos separate contrasts';
write_neg_pos.labels = {'Yes', 'No'};
write_neg_pos.values = {1,0};
write_neg_pos.val    = {0};
write_neg_pos.help = {'If generating negative contrasts, whether to output '
    'separate maps for negative and positive contrasts and for both, '
    'or only the maps with both contrasts' }';


% Executable Branch
liom_contrast      = cfg_exbranch;      
liom_contrast.name = 'Liom Contrast Calculations';            
liom_contrast.tag  = 'liom_contrast';
%PP removed: liom_contrast_struct
liom_contrast.val  = {NIRSmat NewDirCopyNIRS ProcessContrastsBySession GroupMultiSession view consess ...
    spatial_LPF GenerateInverted GroupColorbars contrast_p_value ...
    contrast_figures override_colorbar figures_visible GroupFiguresIntoSubplots ...
    output_unc SmallFigures write_neg_pos TopoData}; %Study_type
liom_contrast.prog = @nirs_run_liom_contrast;  
liom_contrast.vout = @nirs_cfg_vout_liom_contrast; 
liom_contrast.help = {'Liom Contrast Calculations.'};

function vout = nirs_cfg_vout_liom_contrast(job)
vout = cfg_dep;                     
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % NIRS_SPM Group Level Model Estimation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Contrast_files         = cfg_files;  
% Contrast_files.name    = 'Select estimated contrasts files'; 
% Contrast_files.tag     = 'Contrast_files';       
% Contrast_files.filter  = 'mat';    
% Contrast_files.num     = [1 Inf];     
% Contrast_files.help    = {'Select estimated constrast files for this '
%             'group. Select all desired files, and code will try to group '
%             'them by view type, contrast type, and by chromophore. '
%             'Please see code if any doubt.'}'; 
% 
% % Executable Branch
% NIRS_SPM_group      = cfg_exbranch;       
% NIRS_SPM_group.name = 'NIRS_SPM Group Model Estimation';             
% NIRS_SPM_group.tag  = 'NIRS_SPM_group'; 
% NIRS_SPM_group.val  = {Contrast_files}; 
% NIRS_SPM_group.prog = @nirs_run_NIRS_SPM_group;  
% NIRS_SPM_group.vout = @nirs_cfg_vout_NIRS_SPM_group; 
% NIRS_SPM_group.help = {'NIRS_SPM Group level model estimation.'};
% 
% function vout = nirs_cfg_vout_NIRS_SPM_group(job)
% vout = cfg_dep;                     
% vout.sname      = 'NIRS.mat';       
% vout.src_output = substruct('.','NIRSmat'); 
% vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
% end

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



%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPM factorial design configuration
%%%%%%%%%%%%%%%%

% ---------------------------------------------------------------------
% dir Directory
% ---------------------------------------------------------------------
dir         = cfg_files;
dir.tag     = 'dir';
dir.name    = 'Directory';
dir.help    = {'Select a directory where the SPM.mat file containing the specified design matrix will be written.'};
dir.filter = 'dir';
dir.val{1} = {''};
dir.ufilter = '.*';
dir.num     = [0 1];

% ---------------------------------------------------------------------
% scans Scans
% ---------------------------------------------------------------------
scans         = cfg_files;
scans.tag     = 'scans';
scans.name    = 'Scans';
scans.help    = {'Select the images.  They must all have the same image dimensions, orientation, voxel size etc.'};
scans.filter = 'image';
scans.ufilter = '.*';
scans.num     = [0 Inf];
% ---------------------------------------------------------------------
% t1 One-sample t-test
% ---------------------------------------------------------------------
t1         = cfg_branch;
t1.tag     = 't1';
t1.name    = 'One-sample t-test';
t1.val     = {scans };
t1.help    = {''};

% ---------------------------------------------------------------------
% scans1 Group 1 scans
% ---------------------------------------------------------------------
scans1         = cfg_files;
scans1.tag     = 'scans1';
scans1.name    = 'Group 1 scans';
scans1.help    = {'Select the images from sample 1.  They must all have the same image dimensions, orientation, voxel size etc.'};
scans1.filter = 'image';
scans1.ufilter = '.*';
scans1.num     = [0 Inf];
% ---------------------------------------------------------------------
% scans2 Group 2 scans
% ---------------------------------------------------------------------
scans2         = cfg_files;
scans2.tag     = 'scans2';
scans2.name    = 'Group 2 scans';
scans2.help    = {'Select the images from sample 2.  They must all have the same image dimensions, orientation, voxel size etc.'};
scans2.filter = 'image';
scans2.ufilter = '.*';
scans2.num     = [0 Inf];
% ---------------------------------------------------------------------
% dept Independence
% ---------------------------------------------------------------------
dept         = cfg_menu;
dept.tag     = 'dept';
dept.name    = 'Independence';
dept.help    = {
                'By default, the measurements are assumed to be independent between levels. '
                ''
                'If you change this option to allow for dependencies, this will violate the assumption of sphericity. It would therefore be an example of non-sphericity. One such example would be where you had repeated measurements from the same subjects - it may then be the case that, over subjects, measure 1 is correlated to measure 2. '
                ''
                'Restricted Maximum Likelihood (REML): The ensuing covariance components will be estimated using ReML in spm_spm (assuming the same for all responsive voxels) and used to adjust the statistics and degrees of freedom during inference. By default spm_spm will use weighted least squares to produce Gauss-Markov or Maximum likelihood estimators using the non-sphericity structure specified at this stage. The components will be found in SPM.xVi and enter the estimation procedure exactly as the serial correlations in fMRI models.'
                ''
}';
dept.labels  = {
               'Yes'
               'No'
}';
dept.values  = {0 1};
dept.val     = {0};
% ---------------------------------------------------------------------
% deptn Independence (default is 'No')
% ---------------------------------------------------------------------
deptn         = cfg_menu;
deptn.tag     = 'dept';
deptn.name    = 'Independence';
deptn.help    = {
                'By default, the measurements are assumed to be dependent between levels. '
                ''
                'If you change this option to allow for dependencies, this will violate the assumption of sphericity. It would therefore be an example of non-sphericity. One such example would be where you had repeated measurements from the same subjects - it may then be the case that, over subjects, measure 1 is correlated to measure 2. '
                ''
                'Restricted Maximum Likelihood (REML): The ensuing covariance components will be estimated using ReML in spm_spm (assuming the same for all responsive voxels) and used to adjust the statistics and degrees of freedom during inference. By default spm_spm will use weighted least squares to produce Gauss-Markov or Maximum likelihood estimators using the non-sphericity structure specified at this stage. The components will be found in SPM.xVi and enter the estimation procedure exactly as the serial correlations in fMRI models.'
                ''
}';
deptn.labels  = {
               'Yes'
               'No'
}';
deptn.values  = {0 1};
deptn.val     = {1};

% ---------------------------------------------------------------------
% variance Variance
% ---------------------------------------------------------------------
variance         = cfg_menu;
variance.tag     = 'variance';
variance.name    = 'Variance';
variance.help    = {
                    'By default, the measurements in each level are assumed to have unequal variance. '
                    ''
                    'This violates the assumption of ''sphericity'' and is therefore an example of ''non-sphericity''.'
                    ''
                    'This can occur, for example, in a 2nd-level analysis of variance, one contrast may be scaled differently from another.  Another example would be the comparison of qualitatively different dependent variables (e.g. normals vs. patients).  Different variances (heteroscedasticy) induce different error covariance components that are estimated using restricted maximum likelihood (see below).'
                    ''
                    'Restricted Maximum Likelihood (REML): The ensuing covariance components will be estimated using ReML in spm_spm (assuming the same for all responsive voxels) and used to adjust the statistics and degrees of freedom during inference. By default spm_spm will use weighted least squares to produce Gauss-Markov or Maximum likelihood estimators using the non-sphericity structure specified at this stage. The components will be found in SPM.xVi and enter the estimation procedure exactly as the serial correlations in fMRI models.'
                    ''
}';
variance.labels = {
                   'Equal'
                   'Unequal'
}';
variance.values = {0 1};
variance.val    = {1};
% ---------------------------------------------------------------------
% gmsca Grand mean scaling
% ---------------------------------------------------------------------
gmsca         = cfg_menu;
gmsca.tag     = 'gmsca';
gmsca.name    = 'Grand mean scaling';
gmsca.help    = {
                 'This option is only used for PET data.'
                 ''
                 'Selecting YES will specify ''grand mean scaling by factor'' which could be eg. ''grand mean scaling by subject'' if the factor is ''subject''. '
                 ''
                 'Since differences between subjects may be due to gain and sensitivity effects, AnCova by subject could be combined with "grand mean scaling by subject" to obtain a combination of between subject proportional scaling and within subject AnCova. '
                 ''
}';
gmsca.labels = {
                'No'
                'Yes'
}';
gmsca.values = {0 1};
gmsca.val    = {0};
% ---------------------------------------------------------------------
% ancova ANCOVA
% ---------------------------------------------------------------------
ancova         = cfg_menu;
ancova.tag     = 'ancova';
ancova.name    = 'ANCOVA';
ancova.help    = {
                  'This option is only used for PET data.'
                  ''
                  'Selecting YES will specify ''ANCOVA-by-factor'' regressors. This includes eg. ''Ancova by subject'' or ''Ancova by effect''. These options allow eg. different subjects to have different relationships between local and global measurements. '
                  ''
}';
ancova.labels = {
                 'No'
                 'Yes'
}';
ancova.values = {0 1};
ancova.val    = {0};
% ---------------------------------------------------------------------
% t2 Two-sample t-test
% ---------------------------------------------------------------------
t2         = cfg_branch;
t2.tag     = 't2';
t2.name    = 'Two-sample t-test';
t2.val     = {scans1 scans2 dept variance gmsca ancova };
t2.help    = {''};

% ---------------------------------------------------------------------
% scans Scans [1,2]
% ---------------------------------------------------------------------
scans         = cfg_files;
scans.tag     = 'scans';
scans.name    = 'Scans [1,2]';
scans.help    = {'Select the pair of images. '};
scans.filter = 'image';
scans.ufilter = '.*';
scans.num     = [2 2];
% ---------------------------------------------------------------------
% pair Pair
% ---------------------------------------------------------------------
pair         = cfg_branch;
pair.tag     = 'pair';
pair.name    = 'Pair';
pair.val     = {scans };
pair.help    = {'Add a new pair of scans to your experimental design'};
% ---------------------------------------------------------------------
% generic Pairs
% ---------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Pairs';
generic.help    = {''};
generic.values  = {pair};
generic.num     = [1 Inf];
% ---------------------------------------------------------------------
% pt Paired t-test
% ---------------------------------------------------------------------
pt         = cfg_branch;
pt.tag     = 'pt';
pt.name    = 'Paired t-test';
pt.val     = {generic gmsca ancova};
pt.help    = {''};

% ---------------------------------------------------------------------
% scans Scans
% ---------------------------------------------------------------------
scans         = cfg_files;
scans.tag     = 'scans';
scans.name    = 'Scans';
scans.help    = {'Select the images.  They must all have the same image dimensions, orientation, voxel size etc.'};
scans.filter = 'image';
scans.ufilter = '.*';
scans.num     = [1 Inf];
% ---------------------------------------------------------------------
% c Vector
% ---------------------------------------------------------------------
c         = cfg_entry;
c.tag     = 'c';
c.name    = 'Vector';
c.help    = {'Vector of covariate values'};
c.strtype = 'e';
c.num     = [Inf 1];
% ---------------------------------------------------------------------
% cname Name
% ---------------------------------------------------------------------
cname         = cfg_entry;
cname.tag     = 'cname';
cname.name    = 'Name';
cname.help    = {'Name of covariate'};
cname.strtype = 's';
cname.num     = [1 Inf];
% ---------------------------------------------------------------------
% iCC Centering
% ---------------------------------------------------------------------
iCC         = cfg_menu;
iCC.tag     = 'iCC';
iCC.name    = 'Centering';
iCC.help    = {''};
iCC.labels = {
              'Overall mean'
              'No centering'
}';
iCC.values = {1 5};
iCC.val    = {1};
% ---------------------------------------------------------------------
% mcov Covariate
% ---------------------------------------------------------------------
mcov         = cfg_branch;
mcov.tag     = 'mcov';
mcov.name    = 'Covariate';
mcov.val     = {c cname iCC };
mcov.help    = {'Add a new covariate to your experimental design'};
% ---------------------------------------------------------------------
% generic Covariates
% ---------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Covariates';
generic.help    = {'Covariates'};
generic.values  = {mcov };
generic.num     = [0 Inf];
% ---------------------------------------------------------------------
% incint Intercept
% ---------------------------------------------------------------------
incint = cfg_menu;
incint.tag = 'incint';
incint.name = 'Intercept';
incint.help = {['By default, an intercept is always added to the model. If the ',...
    'covariates supplied by the user include a constant effect, the ',...
    'intercept may be omitted.']};
incint.labels = {'Include Intercept','Omit Intercept'};
incint.values = {1,0};
incint.val    = {1};
% ---------------------------------------------------------------------
% mreg Multiple regression
% ---------------------------------------------------------------------
mreg         = cfg_branch;
mreg.tag     = 'mreg';
mreg.name    = 'Multiple regression';
mreg.val     = {scans generic incint};
mreg.help    = {''};
% ---------------------------------------------------------------------
% name Name
% ---------------------------------------------------------------------
name         = cfg_entry;
name.tag     = 'name';
name.name    = 'Name';
name.help    = {'Name of factor, eg. ''Repetition'' '};
name.strtype = 's';
name.num     = [1 Inf];
% ---------------------------------------------------------------------
% levels Levels
% ---------------------------------------------------------------------
levels         = cfg_entry;
levels.tag     = 'levels';
levels.name    = 'Levels';
levels.help    = {'Enter number of levels for this factor, eg. 2'};
levels.strtype = 'e';
levels.num     = [Inf 1];
% ---------------------------------------------------------------------
% fact Factor
% ---------------------------------------------------------------------
fact         = cfg_branch;
fact.tag     = 'fact';
fact.name    = 'Factor';
fact.val     = {name levels dept variance gmsca ancova };
fact.help    = {'Add a new factor to your experimental design'};
% ---------------------------------------------------------------------
% generic Factors
% ---------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Factors';
generic.help    = {
                   'Specify your design a factor at a time. '
                   ''
}';
generic.values  = {fact };
generic.num     = [1 Inf];
% ---------------------------------------------------------------------
% levels Levels
% ---------------------------------------------------------------------
levels         = cfg_entry;
levels.tag     = 'levels';
levels.name    = 'Levels';
levels.help    = {
                  'Enter a vector or scalar that specifies which cell in the factorial design these images belong to. The length of this vector should correspond to the number of factors in the design'
                  ''
                  'For example, length 2 vectors should be used for two-factor designs eg. the vector [2 3] specifies the cell corresponding to the 2nd-level of the first factor and the 3rd level of the 2nd factor.'
                  ''
}';
levels.strtype = 'e';
levels.num     = [Inf 1];
% ---------------------------------------------------------------------
% scans Scans
% ---------------------------------------------------------------------
scans         = cfg_files;
scans.tag     = 'scans';
scans.name    = 'Scans';
scans.help    = {'Select the images for this cell.  They must all have the same image dimensions, orientation, voxel size etc.'};
scans.filter = 'image';
scans.ufilter = '.*';
scans.num     = [1 Inf];
% ---------------------------------------------------------------------
% icell Cell
% ---------------------------------------------------------------------
icell         = cfg_branch;
icell.tag     = 'icell';
icell.name    = 'Cell';
icell.val     = {levels scans };
icell.help    = {'Enter data for a cell in your design'};
% ---------------------------------------------------------------------
% scell Cell
% ---------------------------------------------------------------------
scell         = cfg_branch;
scell.tag     = 'icell';
scell.name    = 'Cell';
scell.val     = {scans };
scell.help    = {'Enter data for a cell in your design'};
% ---------------------------------------------------------------------
% generic Specify cells
% ---------------------------------------------------------------------
generic1         = cfg_repeat;
generic1.tag     = 'generic';
generic1.name    = 'Specify cells';
generic1.help    = {
                    'Enter the scans a cell at a time'
                    ''
}';
generic1.values  = {icell };
generic1.num     = [1 Inf];
% ---------------------------------------------------------------------
% generic Specify cells
% ---------------------------------------------------------------------
generic2         = cfg_repeat;
generic2.tag     = 'generic';
generic2.name    = 'Specify cells';
generic2.help    = {
                    'Enter the scans a cell at a time'
                    ''
}';
generic2.values  = {scell };
generic2.num     = [1 Inf];
% ---------------------------------------------------------------------
% anova ANOVA 
% ---------------------------------------------------------------------
anova         = cfg_branch;
anova.tag     = 'anova';
anova.name    = 'One-way ANOVA';
anova.val     = {generic2 dept variance gmsca ancova};
anova.help    = {
              'One-way Analysis of Variance (ANOVA)'
}';
% ---------------------------------------------------------------------
% fd Full factorial
% ---------------------------------------------------------------------
fd         = cfg_branch;
fd.tag     = 'fd';
fd.name    = 'Full factorial';
fd.val     = {generic generic1 };
fd.help    = {
              'This option is best used when you wish to test for all main effects and interactions in one-way, two-way or three-way ANOVAs. Design specification proceeds in 2 stages. Firstly, by creating new factors and specifying the number of levels and name for each. Nonsphericity, ANOVA-by-factor and scaling options can also be specified at this stage. Secondly, scans are assigned separately to each cell. This accomodates unbalanced designs.'
              ''
              'For example, if you wish to test for a main effect in the population from which your subjects are drawn and have modelled that effect at the first level using K basis functions (eg. K=3 informed basis functions) you can use a one-way ANOVA with K-levels. Create a single factor with K levels and then assign the data to each cell eg. canonical, temporal derivative and dispersion derivative cells, where each cell is assigned scans from multiple subjects.'
              ''
              'SPM will also automatically generate the contrasts necessary to test for all main effects and interactions. '
              ''
}';
% ---------------------------------------------------------------------
% name Name
% ---------------------------------------------------------------------
name         = cfg_entry;
name.tag     = 'name';
name.name    = 'Name';
name.help    = {'Name of factor, eg. ''Repetition'' '};
name.strtype = 's';
name.num     = [1 Inf];
% ---------------------------------------------------------------------
% fac Factor
% ---------------------------------------------------------------------
fac         = cfg_branch;
fac.tag     = 'fac';
fac.name    = 'Factor';
fac.val     = {name dept variance gmsca ancova };
fac.help    = {
               'Add a new factor to your design.'
               ''
               'If you are using the ''Subjects'' option to specify your scans and conditions, you may wish to make use of the following facility. There are two reserved words for the names of factors. These are ''subject'' and ''repl'' (standing for replication). If you use these factor names then SPM can automatically create replication and/or subject factors without you having to type in an extra entry in the condition vector.'
               ''
               'For example, if you wish to model Subject and Task effects (two factors), under Subjects->Subject->Conditions you can type in simply [1 2 1 2] to specify eg. just the ''Task'' factor level. You do not need to eg. for the 4th subject enter the matrix [1 4; 2 4; 1 4; 2 4]. '
               ''
}';
% ---------------------------------------------------------------------
% generic Factors
% ---------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Factors';
generic.help    = {
                   'Specify your design a factor at a time.'
                   ''
}';
generic.values  = {fac };
generic.num     = [1 Inf];
% ---------------------------------------------------------------------
% scans Scans
% ---------------------------------------------------------------------
scans         = cfg_files;
scans.tag     = 'scans';
scans.name    = 'Scans';
scans.help    = {'Select the images to be analysed.  They must all have the same image dimensions, orientation, voxel size etc.'};
scans.filter = 'image';
scans.ufilter = '.*';
scans.num     = [1 Inf];
% ---------------------------------------------------------------------
% conds Conditions
% ---------------------------------------------------------------------
conds         = cfg_entry;
conds.tag     = 'conds';
conds.name    = 'Conditions';
conds.help    = {''};
conds.strtype = 'e';
conds.num     = [Inf Inf];
% ---------------------------------------------------------------------
% fsubject Subject
% ---------------------------------------------------------------------
fsubject         = cfg_branch;
fsubject.tag     = 'fsubject';
fsubject.name    = 'Subject';
fsubject.val     = {scans conds };
fsubject.help    = {'Enter data and conditions for a new subject'};
% ---------------------------------------------------------------------
% generic Subjects
% ---------------------------------------------------------------------
generic1         = cfg_repeat;
generic1.tag     = 'generic';
generic1.name    = 'Subjects';
generic1.help    = {''};
generic1.values  = {fsubject };
generic1.num     = [1 Inf];
% ---------------------------------------------------------------------
% scans Scans
% ---------------------------------------------------------------------
scans         = cfg_files;
scans.tag     = 'scans';
scans.name    = 'Scans';
scans.help    = {'Select the images to be analysed.  They must all have the same image dimensions, orientation, voxel size etc.'};
scans.filter = 'image';
scans.ufilter = '.*';
scans.num     = [1 Inf];
% ---------------------------------------------------------------------
% imatrix Factor matrix
% ---------------------------------------------------------------------
imatrix         = cfg_entry;
imatrix.tag     = 'imatrix';
imatrix.name    = 'Factor matrix';
imatrix.help    = {'Specify factor/level matrix as a nscan-by-4 matrix. Note that the first column of I is reserved for the internal replication factor and must not be used for experimental factors.'};
imatrix.strtype = 'e';
imatrix.num     = [Inf Inf];
% ---------------------------------------------------------------------
% specall Specify all
% ---------------------------------------------------------------------
specall         = cfg_branch;
specall.tag     = 'specall';
specall.name    = 'Specify all';
specall.val     = {scans imatrix };
specall.help    = {
                   'Specify (i) all scans in one go and (ii) all conditions using a factor matrix, I. This option is for ''power users''. The matrix I must have four columns and as as many rows as scans. It has the same format as SPM''s internal variable SPM.xX.I. '
                   ''
                   'The first column of I denotes the replication number and entries in the other columns denote the levels of each experimental factor.'
                   ''
                   'So, for eg. a two-factor design the first column denotes the replication number and columns two and three have entries like 2 3 denoting the 2nd level of the first factor and 3rd level of the second factor. The 4th column in I would contain all 1s.'
}';
% ---------------------------------------------------------------------
% fsuball Specify Subjects or all Scans & Factors
% ---------------------------------------------------------------------
fsuball         = cfg_choice;
fsuball.tag     = 'fsuball';
fsuball.name    = 'Specify Subjects or all Scans & Factors';
fsuball.val     = {generic1 };
fsuball.help    = {''};
fsuball.values  = {generic1 specall };
% ---------------------------------------------------------------------
% fnum Factor number
% ---------------------------------------------------------------------
fnum         = cfg_entry;
fnum.tag     = 'fnum';
fnum.name    = 'Factor number';
fnum.help    = {'Enter the number of the factor.'};
fnum.strtype = 'e';
fnum.num     = [1 1];
% ---------------------------------------------------------------------
% fmain Main effect
% ---------------------------------------------------------------------
fmain         = cfg_branch;
fmain.tag     = 'fmain';
fmain.name    = 'Main effect';
fmain.val     = {fnum };
fmain.help    = {'Add a main effect to your design matrix'};
% ---------------------------------------------------------------------
% fnums Factor numbers
% ---------------------------------------------------------------------
fnums         = cfg_entry;
fnums.tag     = 'fnums';
fnums.name    = 'Factor numbers';
fnums.help    = {'Enter the numbers of the factors of this (two-way) interaction.'};
fnums.strtype = 'e';
fnums.num     = [2 1];
% ---------------------------------------------------------------------
% inter Interaction
% ---------------------------------------------------------------------
inter         = cfg_branch;
inter.tag     = 'inter';
inter.name    = 'Interaction';
inter.val     = {fnums };
inter.help    = {'Add an interaction to your design matrix'};
% ---------------------------------------------------------------------
% maininters Main effects & Interactions
% ---------------------------------------------------------------------
maininters         = cfg_repeat;
maininters.tag     = 'maininters';
maininters.name    = 'Main effects & Interactions';
maininters.help    = {''};
maininters.values  = {fmain inter };
maininters.num     = [1 Inf];
% ---------------------------------------------------------------------
% anovaw ANOVA within subject
% ---------------------------------------------------------------------
anovaw         = cfg_branch;
anovaw.tag     = 'anovaw';
anovaw.name    = 'One-way ANOVA - within subject';
anovaw.val     = {generic1 deptn variance gmsca ancova};
anovaw.help    = {
              'One-way Analysis of Variance (ANOVA) - within subject'
}';
% ---------------------------------------------------------------------
% fblock Flexible factorial
% ---------------------------------------------------------------------
fblock         = cfg_branch;
fblock.tag     = 'fblock';
fblock.name    = 'Flexible factorial';
fblock.val     = {generic fsuball maininters };
fblock.help    = {
                  'Create a design matrix a block at a time by specifying which main effects and interactions you wish to be included.'
                  ''
                  'This option is best used for one-way, two-way or three-way ANOVAs but where you do not wish to test for all possible main effects and interactions. This is perhaps most useful for PET where there is usually not enough data to test for all possible effects. Or for 3-way ANOVAs where you do not wish to test for all of the two-way interactions. A typical example here would be a group-by-drug-by-task analysis where, perhaps, only (i) group-by-drug or (ii) group-by-task interactions are of interest. In this case it is only necessary to have two-blocks in the design matrix - one for each interaction. The three-way interaction can then be tested for using a contrast that computes the difference between (i) and (ii).'
                  ''
                  'Design specification then proceeds in 3 stages. Firstly, factors are created and names specified for each. Nonsphericity, ANOVA-by-factor and scaling options can also be specified at this stage.'
                  ''
                  'Secondly, a list of scans is produced along with a factor matrix, I. This is an nscan x 4 matrix of factor level indicators (see xX.I below). The first factor must be ''replication'' but the other factors can be anything. Specification of I and the scan list can be achieved in one of two ways (a) the ''Specify All'' option allows I to be typed in at the user interface or (more likely) loaded in from the matlab workspace. All of the scans are then selected in one go. (b) the ''Subjects'' option allows you to enter scans a subject at a time. The corresponding experimental conditions (ie. levels of factors) are entered at the same time. SPM will then create the factor matrix I. This style of interface is similar to that available in SPM2.'
                  ''
                  'Thirdly, the design matrix is built up a block at a time. Each block can be a main effect or a (two-way) interaction. '
                  ''
}';
% ---------------------------------------------------------------------
% des Design
% ---------------------------------------------------------------------
des         = cfg_choice;
des.tag     = 'des';
des.name    = 'Design';
des.val     = {t1 };
des.help    = {''};
des.values  = {t1 t2 pt mreg anova anovaw fd fblock };
% ---------------------------------------------------------------------
% c Vector
% ---------------------------------------------------------------------
c         = cfg_entry;
c.tag     = 'c';
c.name    = 'Vector';
c.help    = {
             'Vector of covariate values.'
             'Enter the covariate values ''''per subject'''' (i.e. all for subject 1, then all for subject 2, etc). Importantly, the ordering of the cells of a factorial design has to be the same for all subjects in order to be consistent with the ordering of the covariate values.'
}';
c.strtype = 'e';
c.num     = [Inf 1];
% ---------------------------------------------------------------------
% cname Name
% ---------------------------------------------------------------------
cname         = cfg_entry;
cname.tag     = 'cname';
cname.name    = 'Name';
cname.help    = {'Name of covariate'};
cname.strtype = 's';
cname.num     = [1 Inf];
% ---------------------------------------------------------------------
% iCFI Interactions
% ---------------------------------------------------------------------
iCFI         = cfg_menu;
iCFI.tag     = 'iCFI';
iCFI.name    = 'Interactions';
iCFI.help    = {
                'For each covariate you have defined, there is an opportunity to create an additional regressor that is the interaction between the covariate and a chosen experimental factor. '
                ''
}';
iCFI.labels = {
               'None'
               'With Factor 1'
               'With Factor 2'
               'With Factor 3'
}';
iCFI.values = {1 2 3 4};
iCFI.val    = {1};
% ---------------------------------------------------------------------
% iCC Centering
% ---------------------------------------------------------------------
iCC         = cfg_menu;
iCC.tag     = 'iCC';
iCC.name    = 'Centering';
iCC.help    = {
               'The appropriate centering option is usually the one that corresponds to the interaction chosen, and ensures that main effects of the interacting factor aren''t affected by the covariate. You are advised to choose this option, unless you have other modelling considerations. '
               ''
}';
iCC.labels = {
              'Overall mean'
              'Factor 1 mean'
              'Factor 2 mean'
              'Factor 3 mean'
              'No centering'
              'User specified value'
              'As implied by ANCOVA'
              'GM'
}';
iCC.values = {1 2 3 4 5 6 7 8};
iCC.val    = {1};
% ---------------------------------------------------------------------
% cov Covariate
% ---------------------------------------------------------------------
cov         = cfg_branch;
cov.tag     = 'cov';
cov.name    = 'Covariate';
cov.val     = {c cname iCFI iCC };
cov.help    = {'Add a new covariate to your experimental design'};
% ---------------------------------------------------------------------
% generic Covariates
% ---------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Covariates';
generic.help    = {
                   'This option allows for the specification of covariates and nuisance variables. Unlike SPM94/5/6, where the design was partitioned into effects of interest and nuisance effects for the computation of adjusted data and the F-statistic (which was used to thresh out voxels where there appeared to be no effects of interest), SPM does not partition the design in this way anymore. The only remaining distinction between effects of interest (including covariates) and nuisance effects is their location in the design matrix, which we have retained for continuity.  Pre-specified design matrix partitions can be entered. '
                   ''
}';
generic.values  = {cov };
generic.num     = [0 Inf];
% ---------------------------------------------------------------------
% tm_none None
% ---------------------------------------------------------------------
tm_none         = cfg_const;
tm_none.tag     = 'tm_none';
tm_none.name    = 'None';
tm_none.val     = {1};
tm_none.help    = {'No threshold masking'};
% ---------------------------------------------------------------------
% athresh Threshold
% ---------------------------------------------------------------------
athresh         = cfg_entry;
athresh.tag     = 'athresh';
athresh.name    = 'Threshold';
athresh.help    = {
                   'Enter the absolute value of the threshold.'
                   ''
}';
athresh.strtype = 'e';
athresh.num     = [1 1];
athresh.val     = {100};
% ---------------------------------------------------------------------
% tma Absolute
% ---------------------------------------------------------------------
tma         = cfg_branch;
tma.tag     = 'tma';
tma.name    = 'Absolute';
tma.val     = {athresh };
tma.help    = {
               'Images are thresholded at a given value and only voxels at which all images exceed the threshold are included. '
               ''
               'This option allows you to specify the absolute value of the threshold.'
               ''
}';
% ---------------------------------------------------------------------
% rthresh Threshold
% ---------------------------------------------------------------------
rthresh         = cfg_entry;
rthresh.tag     = 'rthresh';
rthresh.name    = 'Threshold';
rthresh.help    = {
                   'Enter the threshold as a proportion of the global value'
                   ''
}';
rthresh.strtype = 'e';
rthresh.num     = [1 1];
rthresh.val     = {.8};
% ---------------------------------------------------------------------
% tmr Relative
% ---------------------------------------------------------------------
tmr         = cfg_branch;
tmr.tag     = 'tmr';
tmr.name    = 'Relative';
tmr.val     = {rthresh };
tmr.help    = {
               'Images are thresholded at a given value and only voxels at which all images exceed the threshold are included. '
               ''
               'This option allows you to specify the value of the threshold as a proportion of the global value. '
               ''
}';
% ---------------------------------------------------------------------
% tm Threshold masking
% ---------------------------------------------------------------------
tm         = cfg_choice;
tm.tag     = 'tm';
tm.name    = 'Threshold masking';
tm.val     = {tm_none };
tm.help    = {
              'Images are thresholded at a given value and only voxels at which all images exceed the threshold are included. '
              ''
}';
tm.values  = {tm_none tma tmr };
% ---------------------------------------------------------------------
% im Implicit Mask
% ---------------------------------------------------------------------
im         = cfg_menu;
im.tag     = 'im';
im.name    = 'Implicit Mask';
im.help    = {
              'An "implicit mask" is a mask implied by a particular voxel value. Voxels with this mask value are excluded from the analysis. '
              ''
              'For image data-types with a representation of NaN (see spm_type.m), NaN''s is the implicit mask value, (and NaN''s are always masked out). '
              ''
              'For image data-types without a representation of NaN, zero is the mask value, and the user can choose whether zero voxels should be masked out or not.'
              ''
              'By default, an implicit mask is used. '
              ''
}';
im.labels = {
             'Yes'
             'No'
}';
im.values = {1 0};
im.val    = {1};
% ---------------------------------------------------------------------
% em Explicit Mask
% ---------------------------------------------------------------------
em         = cfg_files;
em.tag     = 'em';
em.name    = 'Explicit Mask';
em.val     = {{''}};
em.help    = {
              'Explicit masks are other images containing (implicit) masks that are to be applied to the current analysis.'
              ''
              'All voxels with value NaN (for image data-types with a representation of NaN), or zero (for other data types) are excluded from the analysis. '
              ''
              'Explicit mask images can have any orientation and voxel/image size. Nearest neighbour interpolation of a mask image is used if the voxel centers of the input images do not coincide with that of the mask image.'
              ''
}';
em.filter = 'image';
em.ufilter = '.*';
em.num     = [0 1];
% ---------------------------------------------------------------------
% masking Masking
% ---------------------------------------------------------------------
masking         = cfg_branch;
masking.tag     = 'masking';
masking.name    = 'Masking';
masking.val     = {tm im em };
masking.help    = {
                   'The mask specifies the voxels within the image volume which are to be assessed. SPM supports three methods of masking (1) Threshold, (2) Implicit and (3) Explicit. The volume analysed is the intersection of all masks.'
                   ''
}';
% ---------------------------------------------------------------------
% g_omit Omit
% ---------------------------------------------------------------------
g_omit         = cfg_const;
g_omit.tag     = 'g_omit';
g_omit.name    = 'Omit';
g_omit.val     = {1};
g_omit.help    = {'Omit'};
% ---------------------------------------------------------------------
% global_uval Global values
% ---------------------------------------------------------------------
global_uval         = cfg_entry;
global_uval.tag     = 'global_uval';
global_uval.name    = 'Global values';
global_uval.help    = {
                       'Enter the vector of global values'
                       ''
}';
global_uval.strtype = 'e';
global_uval.num     = [Inf 1];
% ---------------------------------------------------------------------
% g_user User
% ---------------------------------------------------------------------
g_user         = cfg_branch;
g_user.tag     = 'g_user';
g_user.name    = 'User';
g_user.val     = {global_uval };
g_user.help    = {
                  'User defined  global effects (enter your own '
                  'vector of global values)'
}';
% ---------------------------------------------------------------------
% g_mean Mean
% ---------------------------------------------------------------------
g_mean         = cfg_const;
g_mean.tag     = 'g_mean';
g_mean.name    = 'Mean';
g_mean.val     = {1};
g_mean.help    = {
                  'SPM standard mean voxel value'
                  ''
                  'This defines the global mean via a two-step process. Firstly, the overall mean is computed. Voxels with values less than 1/8 of this value are then deemed extra-cranial and get masked out. The mean is then recomputed on the remaining voxels.'
                  ''
}';
% ---------------------------------------------------------------------
% globalc Global calculation
% ---------------------------------------------------------------------
globalc         = cfg_choice;
globalc.tag     = 'globalc';
globalc.name    = 'Global calculation';
globalc.val     = {g_omit };
globalc.help    = {
                   'This option is only used for PET data.'
                   ''
                   'There are three methods for estimating global effects (1) Omit (assumming no other options requiring the global value chosen) (2) User defined (enter your own vector of global values) (3) Mean: SPM standard mean voxel value (within per image fullmean/8 mask) '
                   ''
}';
globalc.values  = {g_omit g_user g_mean };
% ---------------------------------------------------------------------
% gmsca_no No
% ---------------------------------------------------------------------
gmsca_no         = cfg_const;
gmsca_no.tag     = 'gmsca_no';
gmsca_no.name    = 'No';
gmsca_no.val     = {1};
gmsca_no.help    = {'No overall grand mean scaling'};
% ---------------------------------------------------------------------
% gmscv Grand mean scaled value
% ---------------------------------------------------------------------
gmscv         = cfg_entry;
gmscv.tag     = 'gmscv';
gmscv.name    = 'Grand mean scaled value';
gmscv.help    = {
                 'The default value of 50, scales the global flow to a physiologically realistic value of 50ml/dl/min.'
                 ''
}';
gmscv.strtype = 'e';
gmscv.num     = [Inf 1];
gmscv.val     = {50};
% ---------------------------------------------------------------------
% gmsca_yes Yes
% ---------------------------------------------------------------------
gmsca_yes         = cfg_branch;
gmsca_yes.tag     = 'gmsca_yes';
gmsca_yes.name    = 'Yes';
gmsca_yes.val     = {gmscv };
gmsca_yes.help    = {
                     'Scaling of the overall grand mean simply scales all the data by a common factor such that the mean of all the global values is the value specified. For qualitative data, this puts the data into an intuitively accessible scale without altering the statistics. '
                     ''
}';
% ---------------------------------------------------------------------
% gmsca Overall grand mean scaling
% ---------------------------------------------------------------------
gmsca         = cfg_choice;
gmsca.tag     = 'gmsca';
gmsca.name    = 'Overall grand mean scaling';
gmsca.val     = {gmsca_no };
gmsca.help    = {
                 'Scaling of the overall grand mean simply scales all the data by a common factor such that the mean of all the global values is the value specified. For qualitative data, this puts the data into an intuitively accessible scale without altering the statistics. '
                 ''
                 'When proportional scaling global normalisation is used each image is separately scaled such that it''s global value is that specified (in which case the grand mean is also implicitly scaled to that value). So, to proportionally scale each image so that its global value is eg. 20, select <Yes> then type in 20 for the grand mean scaled value.'
                 ''
                 'When using AnCova or no global normalisation, with data from different subjects or sessions, an intermediate situation may be appropriate, and you may be given the option to scale group, session or subject grand means separately. '
                 ''
}';
gmsca.values  = {gmsca_no gmsca_yes };
% ---------------------------------------------------------------------
% glonorm Normalisation
% ---------------------------------------------------------------------
glonorm         = cfg_menu;
glonorm.tag     = 'glonorm';
glonorm.name    = 'Normalisation';
glonorm.help    = {
                   'Global nuisance effects are usually accounted for either by scaling the images so that they all have the same global value (proportional scaling), or by including the global covariate as a nuisance effect in the general linear model (AnCova). Much has been written on which to use, and when. Basically, since proportional scaling also scales the variance term, it is appropriate for situations where the global measurement predominantly reflects gain or sensitivity. Where variance is constant across the range of global values, linear modelling in an AnCova approach has more flexibility, since the model is not restricted to a simple proportional regression. '
                   ''
                   '''Ancova by subject'' or ''Ancova by effect'' options are implemented using the ANCOVA options provided where each experimental factor (eg. subject or effect), is defined. These allow eg. different subjects to have different relationships between local and global measurements. '
                   ''
                   'Since differences between subjects may be due to gain and sensitivity effects, AnCova by subject could be combined with "grand mean scaling by subject" (an option also provided where each experimental factor is originally defined) to obtain a combination of between subject proportional scaling and within subject AnCova. '
                   ''
}';
glonorm.labels = {
                  'None'
                  'Proportional'
                  'ANCOVA'
}';
glonorm.values = {1 2 3};
glonorm.val    = {1};
% ---------------------------------------------------------------------
% globalm Global normalisation
% ---------------------------------------------------------------------
globalm         = cfg_branch;
globalm.tag     = 'globalm';
globalm.name    = 'Global normalisation';
globalm.val     = {gmsca glonorm };
globalm.help    = {
                   'This option is only used for PET data.'
                   ''
                   'Global nuisance effects are usually accounted for either by scaling the images so that they all have the same global value (proportional scaling), or by including the global covariate as a nuisance effect in the general linear model (AnCova). Much has been written on which to use, and when. Basically, since proportional scaling also scales the variance term, it is appropriate for situations where the global measurement predominantly reflects gain or sensitivity. Where variance is constant across the range of global values, linear modelling in an AnCova approach has more flexibility, since the model is not restricted to a simple proportional regression. '
                   ''
                   '''Ancova by subject'' or ''Ancova by effect'' options are implemented using the ANCOVA options provided where each experimental factor (eg. subject or effect), is defined. These allow eg. different subjects to have different relationships between local and global measurements. '
                   ''
                   'Since differences between subjects may be due to gain and sensitivity effects, AnCova by subject could be combined with "grand mean scaling by subject" (an option also provided where each experimental factor is originally defined) to obtain a combination of between subject proportional scaling and within subject AnCova. '
                   ''
}';
% ---------------------------------------------------------------------
% factorial_design Factorial design specification
% ---------------------------------------------------------------------
factorial_design         = cfg_branch;
factorial_design.tag     = 'factorial_design';
factorial_design.name    = 'Factorial design specification';
factorial_design.val     = {dir des generic masking globalc globalm };
factorial_design.help    = {
                            'This interface is used for setting up analyses of PET data. It is also used for ''2nd level'' or ''random effects'' analysis which allow one to make a population inference. First level models can be used to produce appropriate summary data, which can then be used as raw data for a second-level analysis. For example, a simple t-test on contrast images from the first-level turns out to be a random-effects analysis with random subject effects, inferring for the population based on a particular sample of subjects.'
                            ''
                            'This interface configures the design matrix, describing the general linear model, data specification, and other parameters necessary for the statistical analysis. These parameters are saved in a configuration file (SPM.mat), which can then be passed on to spm_spm.m which estimates the design. This is achieved by pressing the ''Estimate'' button. Inference on these estimated parameters is then handled by the SPM results section. '
                            ''
                            'A separate interface handles design configuration for fMRI time series.'
                            ''
                            'Various data and parameters need to be supplied to specify the design (1) the image files, (2) indicators of the corresponding condition/subject/group (2) any covariates, nuisance variables, or design matrix partitions (3) the type of global normalisation (if any) (4) grand mean scaling options (5) thresholds and masks defining the image volume to analyse. The interface supports a comprehensive range of options for all these parameters.'
                            ''
}';

%-------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%
group_session_to_average         = cfg_entry;
group_session_to_average.name    = 'Session to average';
group_session_to_average.tag     = 'group_session_to_average';       
group_session_to_average.strtype = 'r';
group_session_to_average.num     = [1 1];   
group_session_to_average.val     = {1};
group_session_to_average.help    = {'This is only used for multi-session group studies.'
    'Specify here which session the group analysis should be done upon. '
    'Only one session can be specified.'}'; 

% Executable Branch
liom_group      = cfg_exbranch;       
liom_group.name = 'Liom Group Model Estimation';             
liom_group.tag  = 'liom_group'; 
liom_group.val  = {NIRSmat FFX_or_RFX contrast_figures contrast_p_value ...
        GenerateInverted GroupColorbars override_colorbar figures_visible ...
        GroupFiguresIntoSubplots output_unc SmallFigures write_neg_pos ...
        group_session_to_average}; % factorial_design}; 
liom_group.prog = @nirs_run_liom_group;  
liom_group.vout = @nirs_cfg_vout_liom_group; 
liom_group.help = {'Liom Group level model estimation.'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
ROCiternum.val     = {10};
ROCiternum.strtype = 'r';
ROCiternum.num     = [1 1];     
ROCiternum.help    = {'Number of iterations'}; 

RunGLMorFigures      = cfg_menu;
RunGLMorFigures.tag  = 'RunGLMorFigures';
RunGLMorFigures.name = 'Run GLMs or generate figures';
RunGLMorFigures.labels = {'GLM','Figures','Both'};
RunGLMorFigures.values = {1,2,3};
RunGLMorFigures.val    = {3};
RunGLMorFigures.help = {'Run GLMs and/or generate figures from'
    'Previously run GLMs.'}';

%ROC options                    
ROCnumCh         = cfg_entry; 
ROCnumCh.name    = 'Total Number of channels';
ROCnumCh.tag     = 'ROCnumCh';  
ROCnumCh.val     = {40};
ROCnumCh.strtype = 'r';
ROCnumCh.num     = [1 1];     
ROCnumCh.help    = {'Total Number of channels'}; 

dir_dataSPM         = cfg_entry; 
dir_dataSPM.name    = 'Directory to work from';
dir_dataSPM.tag     = 'dir_dataSPM';  
dir_dataSPM.val{1}  = 'dataSPM';
dir_dataSPM.strtype = 's';
dir_dataSPM.num     = [1 Inf];     
dir_dataSPM.help    = {'Directory to work from'}; 

Volt2         = cfg_entry; 
Volt2.name    = 'Positive or negative t-test for 2nd Volterra';
Volt2.tag     = 'Volt2';  
Volt2.val     = {0};
Volt2.strtype = 'r';
Volt2.num     = [1 Inf];     
Volt2.help    = {'Positive or negative t-test for 2nd Volterra'
    'Enter a matrix of number of jobs by number of subjects'
    'With entries of 0 if negative t-test and 1 if positive t-test'
    'for the 2nd Volterra'
    'Or enter just 0 (and not a matrix) if all t-tests are negative.'}'; 

byIter           = cfg_menu;
byIter.name      = 'Give test result by iteration';
byIter.tag       = 'byIter';
byIter.labels    = {'No' 'Yes'};
byIter.values    = {0,1};
byIter.val       = {0};
byIter.help      = {'Usually, No.'}';

compute_OR           = cfg_menu;
compute_OR.name      = 'Compute HbO and HbR separately';
compute_OR.tag       = 'compute_OR';
compute_OR.labels    = {'No' 'Yes'};
compute_OR.values    = {0,1};
compute_OR.val       = {0};
compute_OR.help      = {'Usually, No.'}';

compute_LU           = cfg_menu;
compute_LU.name      = 'Compute lower and upper bounds';
compute_LU.tag       = 'compute_LU';
compute_LU.labels    = {'No' 'Yes'};
compute_LU.values    = {0,1};
compute_LU.val       = {0};
compute_LU.help      = {'Usually, No.'}';

runFtest           = cfg_menu;
runFtest.name      = 'Run F test';
runFtest.tag       = 'runFtest';
runFtest.labels    = {'No' 'Yes'};
runFtest.values    = {0,1};
runFtest.val       = {0};
runFtest.help      = {'Usually, No.'
                'Careful, can be very slow, since...'}';

% Executable Branch
ROCtest      = cfg_exbranch;       
ROCtest.name = 'ROC Sensitivity and specificity testing';            
ROCtest.tag  = 'ROCtest';
ROCtest.val  = {NIRSmat ROCLoopJob ROCDeleteLarge ROCiternum RunGLMorFigures ...
    ROCnumCh dir_dataSPM Volt2 byIter compute_OR compute_LU runFtest}; 
ROCtest.prog = @nirs_run_ROCtest;  
ROCtest.vout = @nirs_cfg_vout_ROCtest; 
ROCtest.help = {'This module performs a large loop over GLMs'
            'specified with different random seeds. '
            'To use it, user need to first specify in the Matlabbatch'
            'front end a sequence of modules to be run, typically starting'
            'with a module that requires a random seed (such as the '
            'AddTestStimuli module). This sequence of modules is referred '
            'to as a job. The code will run that job repetitively '
            'by incrementing the random seed as many times as specified. '
            'One can later use the separate script ROCfigures_script'
            'to generate a variety of figures: ROC plots and boxplots.'
            'This script will require adapting it to specific requirements.'}';

function vout = nirs_cfg_vout_ROCtest(job)
vout = cfg_dep;                     
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % NIRS_SPM GLM Results Display
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Executable Branch
% NIRS_SPM_model_display      = cfg_exbranch;       
% NIRS_SPM_model_display.name = 'NIRS_SPM Results Display';            
% NIRS_SPM_model_display.tag  = 'NIRS_SPM_model_display';
% NIRS_SPM_model_display.val  = {NIRS_SPM_Coregistration_Channels}; 
% NIRS_SPM_model_display.prog = @nirs_run_NIRS_SPM_model_display;  
% NIRS_SPM_model_display.vout = @nirs_cfg_vout_NIRS_SPM_model_display; 
% NIRS_SPM_model_display.help = {'NIRS_SPM Results Display.'};
% 
% function vout = nirs_cfg_vout_NIRS_SPM_model_display(job)
% vout = cfg_dep;                     
% vout.sname      = 'NIRS.mat';       
% vout.src_output = substruct('.','NIRSmat'); 
% vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % NIRS_SPM Contrasts Display
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %Select map of interest
map_file         = cfg_files;  
map_file.name    = 'Select statistical map'; % The displayed name
map_file.tag     = 'map_file';       
map_file.filter  = 'mat';    
map_file.num     = [1 1];     % Number of inputs required 
map_file.help    = {'Select statistical map of interest (file '
            'containing SPM_nirs (for group) or cinterp_SPM_nirs '
            'structure (for individual) ) - used to select channel '
            'with nearest projected distance to statistical map maximum.'}'; % help text displayed
% 
% %select data file of interest
% data_file         = cfg_files;  
% data_file.name    = 'Select NIRS_SPM data file'; 
% data_file.tag     = 'data_file';       
% data_file.filter  = 'mat';    
% data_file.num     = [1 1];     
% data_file.help    = {'Select data file to calculate contrast estimates '
%                 '(file containing SPM_nirs structure) -- does not have '
%                 'to be related to statistical map.'}'; 
% 
% %Select contrasts
% reg_num         = cfg_entry; %
% reg_num.name    = 'Regressor identification numbers'; 
% reg_num.tag     = 'reg_num';       
% reg_num.strtype = 'r'; 
% reg_num.num     = [1 Inf];     
% reg_num.help    = {'Enter regressor numbers as a Matlab row vector '
%         '(get from the design matrix associated with the data file)'}'; 
% 
% % Executable Branch
% NIRS_SPM_contrast_display      = cfg_exbranch;       
% NIRS_SPM_contrast_display.name = 'NIRS_SPM Contrasts Estimate Display';             
% NIRS_SPM_contrast_display.tag  = 'NIRS_SPM_contrast_display'; 
% %NIRSmat not used currently
% NIRS_SPM_contrast_display.val  = {map_file data_file ...
%             NIRS_SPM_Coregistration_Channels view reg_num};  
% NIRS_SPM_contrast_display.prog = @nirs_run_NIRS_SPM_contrast_display;  
% NIRS_SPM_contrast_display.vout = @nirs_cfg_vout_NIRS_SPM_contrast_display; 
% NIRS_SPM_contrast_display.help = {'NIRS_SPM Contrast Estimates Results Display.'};
% 
% function vout = nirs_cfg_vout_NIRS_SPM_contrast_display(job)
% vout = cfg_dep;                     
% vout.sname      = 'NIRS.mat';       
% vout.src_output = substruct('.','NIRSmat'); 
% vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuration NIRS_SPM diagnostics for protocole and detrending time-series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %select data files of interest
% data_files         = cfg_files;  
% data_files.name    = 'Select NIRS_SPM data files'; 
% data_files.tag     = 'data_files';       
% data_files.filter  = 'mat';    
% data_files.num     = [1 Inf];    
% data_files.help    = {'Select data files for diagnostic (containing '
%                 'SPM_nirs structure).'}'; 
% 
% ch_num         = cfg_entry; 
% ch_num.name    = 'Channel identification numbers'; 
% ch_num.tag     = 'ch_num';       
% ch_num.strtype = 'r';     
% ch_num.num     = [1 Inf];     
% ch_num.help    = {'Enter channel numbers as a Matlab row vector '
%             '(get from the design matrix associated with the data file)'}'; 
% 
% % Executable Branch
% NIRS_SPM_diagnostic      = cfg_exbranch;       
% NIRS_SPM_diagnostic.name = 'NIRS_SPM diagnostics for protocole and detrending time-series';             
% NIRS_SPM_diagnostic.tag  = 'NIRS_SPM_diagnostic'; 
% NIRS_SPM_diagnostic.val  = {data_files ch_num}; 
% NIRS_SPM_diagnostic.prog = @nirs_run_NIRS_SPM_diagnostic; 
% NIRS_SPM_diagnostic.vout = @nirs_cfg_vout_NIRS_SPM_diagnostic; 
% NIRS_SPM_diagnostic.help = {'NIRS_SPM diagnostics for protocole and detrending time-series.'};
% 
% function vout = nirs_cfg_vout_NIRS_SPM_diagnostic(job)
% vout = cfg_dep;                     
% vout.sname      = 'NIRS.mat';       
% vout.src_output = substruct('.','NIRSmat'); 
% vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Extract map data from group and session analyses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

extract_auto_modality   = cfg_menu;
extract_auto_modality.tag  = 'extract_auto_modality';
extract_auto_modality.name = 'Modality map to extract from';
extract_auto_modality.labels = {'HbO','HbR','HbT'};
extract_auto_modality.values = {1,2,3};
extract_auto_modality.val = {2};
extract_auto_modality.help = {'Modality map to extract from.'};

extract_max_HbR           = cfg_branch;
extract_max_HbR.name      = 'extract_max_HbR';
extract_max_HbR.tag       = 'extract_max_HbR';
extract_max_HbR.val       = {extract_auto_modality}; 
extract_max_HbR.help      = {''};

extract_max_all           = cfg_branch;
extract_max_all.name      = 'extract_max_all';
extract_max_all.tag       = 'extract_max_all';
%extract_max_all.val       = {}; 
extract_max_all.help      = {''};

extract_coord_max         = cfg_entry;
extract_coord_max.name    = 'extract_coord_max';
extract_coord_max.tag     = 'extract_coord_max';       
extract_coord_max.strtype = 'r';
extract_coord_max.num     = [1 2];   
extract_coord_max.val     = {[100 200]};
extract_coord_max.help    = {'Enter location of the maximum in pixels. '
    '[0 0] is the upper left corner of the map.'}'; 

extract_coord_min         = cfg_entry;
extract_coord_min.name    = 'extract_coord_min';
extract_coord_min.tag     = 'extract_coord_min';       
extract_coord_min.strtype = 'r';
extract_coord_min.num     = [1 2];   
extract_coord_min.val     = {[150 250]};
extract_coord_min.help    = {'Enter location of the minimum in pixels. '
    '[0 0] is the upper left corner of the map.'}';

extract_coordinates           = cfg_branch;
extract_coordinates.name      = 'extract_coordinates';
extract_coordinates.tag       = 'extract_coordinates';
extract_coordinates.val       = {extract_coord_max extract_coord_min}; 
extract_coordinates.help      = {''};

extract_select_auto_mode           = cfg_choice;
extract_select_auto_mode.name      = 'Automatic Extraction Selection Mode';
extract_select_auto_mode.tag       = 'extract_select_auto_mode';
extract_select_auto_mode.values    = {extract_max_HbR extract_max_all extract_coordinates}; 
extract_select_auto_mode.val       = {extract_max_HbR}; 
extract_select_auto_mode.help      = {'Automatic Extraction Selection Mode.'}'; 

extract_auto_mode           = cfg_branch;
extract_auto_mode.name      = 'extract_auto_mode';
extract_auto_mode.tag       = 'extract_auto_mode';
extract_auto_mode.val       = {extract_select_auto_mode}; 
extract_auto_mode.help      = {''};

extract_manual_modality   = cfg_menu;
extract_manual_modality.tag  = 'extract_manual_modality';
extract_manual_modality.name = 'Modality map to extract from';
extract_manual_modality.labels = {'HbO','HbR','HbT'};
extract_manual_modality.values = {1,2,3};
extract_manual_modality.val = {2};
extract_manual_modality.help = {'Modality map to extract from.'};

extract_manual_mode           = cfg_branch;
extract_manual_mode.name      = 'extract_manual_mode';
extract_manual_mode.tag       = 'extract_manual_mode';
extract_manual_mode.val       = {extract_manual_modality}; 
extract_manual_mode.help      = {''};

extract_select_mode           = cfg_choice;
extract_select_mode.name      = 'Extraction Selection Mode';
extract_select_mode.tag       = 'extract_select_mode';
extract_select_mode.values    = {extract_auto_mode extract_manual_mode}; 
extract_select_mode.val       = {extract_auto_mode}; 
extract_select_mode.help      = {'Extraction Selection Mode.'}'; 

%extract name
extract_struct_name         = cfg_entry; %path
extract_struct_name.name    = 'extract_struct_name';
extract_struct_name.tag     = 'extract_struct_name';       
extract_struct_name.strtype = 's';
extract_struct_name.num     = [1 Inf];     
extract_struct_name.val     = {'ED'}; 
extract_struct_name.help    = {'Name of structure of extracted data to be saved.'}; 

extract_base_contrast         = cfg_entry;
extract_base_contrast.name    = 'Base contrast number';
extract_base_contrast.tag     = 'extract_base_contrast';       
extract_base_contrast.strtype = 'r';
extract_base_contrast.num     = [1 1];   
extract_base_contrast.val     = {1};
extract_base_contrast.help    = {'Enter base contrast number to extract'
    'maps data from. Specify here the contrast based on which the map points'
    'will be selected.'}';

extract_contrast         = cfg_entry;
extract_contrast.name    = 'List of contrast number(s) to extract';
extract_contrast.tag     = 'extract_contrast';       
extract_contrast.strtype = 'r';
extract_contrast.num     = [1 Inf];   
extract_contrast.val     = {1};
extract_contrast.help    = {'Enter contrast number(s) to extract'
    'maps data from. If an array is entered, a loop over all specified'
    'contrasts will be carried out.'}'; 

extract_threshold_val         = cfg_entry;
extract_threshold_val.name    = 'Threshold value';
extract_threshold_val.tag     = 'extract_threshold_val';       
extract_threshold_val.strtype = 'r';
extract_threshold_val.num     = [1 1];   
extract_threshold_val.val     = {3.5};
extract_threshold_val.help    = {'Threshold value of cluster to average.'}'; 

extract_radius_val         = cfg_entry;
extract_radius_val.name    = 'Radius value';
extract_radius_val.tag     = 'extract_radius_val';       
extract_radius_val.strtype = 'r';
extract_radius_val.num     = [1 1];  
extract_radius_val.val     = {3};
extract_radius_val.help    = {'Radius value in pixels for average.'
    '1 pixel corresponds to very approximately 1 mm.'}'; 

extract_threshold           = cfg_branch;
extract_threshold.name      = 'extract_threshold';
extract_threshold.tag       = 'extract_threshold';
extract_threshold.val       = {extract_threshold_val extract_radius_val}; 
extract_threshold.help      = {'Note that here both a threshold and the radius are used as pixel'
    'selection criteria'}';

extract_radius           = cfg_branch;
extract_radius.name      = 'extract_radius';
extract_radius.tag       = 'extract_radius';
extract_radius.val       = {extract_radius_val}; 
extract_radius.help      = {''};

extract_average_mode           = cfg_choice;
extract_average_mode.name      = 'Extraction Averaging Mode';
extract_average_mode.tag       = 'extract_average_mode';
extract_average_mode.values    = {extract_radius extract_threshold}; 
extract_average_mode.val       = {extract_radius}; 
extract_average_mode.help      = {'Extraction Averaging Mode.'}'; 

Volterra_ratio   = cfg_menu;
Volterra_ratio.tag  = 'Volterra_ratio';
Volterra_ratio.name = 'Calculate Volterra ratio';
Volterra_ratio.labels = {'Yes','No'};
Volterra_ratio.values = {1,0};
Volterra_ratio.val = {0};
Volterra_ratio.help = {'Add ratio of 2nd to 1st Volterra to bNmin, bNmax info.'};

% Executable Branch
extract_map_data      = cfg_exbranch;      
extract_map_data.name = 'Extract map data';            
extract_map_data.tag  = 'extract_map_data';
extract_map_data.val  = {NIRSmat view extract_base_contrast ...
    extract_contrast extract_select_mode ...
    extract_average_mode extract_struct_name Volterra_ratio}; 
extract_map_data.prog = @nirs_run_extract_map_data;  
extract_map_data.vout = @nirs_cfg_vout_extract_map_data; 
extract_map_data.help = {'Extract_map_data.'};

function vout = nirs_cfg_vout_extract_map_data(job)
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
coregNIRS.values = {coreg1 coreg2 coreg_manual1 view3d1 resize1};
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
model_reconstruct.values = {tikhonov1 ReMLreconstruct1 testreconstruct1 checkreconstruct1};
model_reconstruct.help   = {'3D Reconstruction of NIRS data.'};

%module 9
model_specify        = cfg_choice; %cfg_repeat;
model_specify.name   = 'GLM Specification';
model_specify.tag    = 'model_specify';
model_specify.values = {wls_bglm_specify}; 
model_specify.help   = {'These modules specify a GLM.'};

%module 10
model_estimate        = cfg_choice; %cfg_repeat; 
model_estimate.name   = 'GLM Estimation';
model_estimate.tag    = 'model_estimate';
model_estimate.values = {wls_bglm_estimate liom_contrast  ...
            liom_group extract_map_data AnalyzeGLM ROCtest}; 
model_estimate.help   = {'These modules estimate a GLM.'};
 

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
    configMC1 runMC1 makesens1 calculatePVE1 model_reconstruct model_specify ...
    model_estimate NIRS_HDM CRIUGM}; %model_display
end