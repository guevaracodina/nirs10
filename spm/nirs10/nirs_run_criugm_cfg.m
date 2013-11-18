function criugm1 = nirs_run_criugm_cfg
redo1      = cfg_menu;
redo1.tag  = 'force_redo';
redo1.name = 'Force processing';
redo1.labels = {'False','True'};
redo1.values = {0,1};
redo1.val  = {0};
redo1.help = {'Force redoing this processing even when it has been done already'};

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

%%%
choice_biopac       = cfg_entry;
choice_biopac.tag   = 'choice_biopac';
choice_biopac.name  = 'Number of ports from biopac module';
choice_biopac.num   = [1 1];
choice_biopac.help  = {'Select number of the auxilary files from Biopac to be loaded.'};

no_biopac      = cfg_branch;
no_biopac.name = 'No Biopac recording';
no_biopac.tag  = 'no_biopac';
no_biopac.help = {'All aux data from ''.nirs'' file. will be considered as TTL from e-prime'};

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
helm_temp.filter  = 'mat';
helm_temp.ufilter = 'NIRS.mat';
helm_temp.val{1}  = {''};
helm_temp.num     = [0 1];
helm_temp.help = {['If you have chosen before ''template'' in choice : ''Individual T1 or template''.'...
    'If you have generated a template for a special helmet and that you want to coregister it with the subject anatomical image, please choose the NIRS.mat you have generated.']};

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
CWsystem.labels = {'CW5','CW6','ImagincV2'};
CWsystem.values = {5,6,2};
CWsystem.def  = @(val)nirs_get_defaults('readNIRS.criugm1.CWsystem', val{:});
CWsystem.help = {'ImagincV2: 16 sources + 16 detectors, no auxiliaries'};

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

study_path        = cfg_files;
study_path.tag     = 'study_path';
study_path.name    = 'Choose study path';
study_path.filter  = 'dir';
study_path.ufilter = '.*';
study_path.num     = [1 1];
study_path.help    = {'Choose directory where you want to put your study. If the directory does not exist you must create it then choose it in the batch.'};

indvdata_chosen      = cfg_branch;
indvdata_chosen.name = 'One set of data per subject';
indvdata_chosen.tag  = 'indvdata_chosen';
indvdata_chosen.help = {'You will have to choose one T1 image in the field ''Raw anatomical image'' and one helmet in the field ''Helmet->Text file from Brainsight.''.'};

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

%%%
biopac         = cfg_choice;
biopac.tag     = 'biopac';
biopac.name    = 'Biopac';
biopac.values  = {choice_biopac no_biopac};
biopac.val     = {no_biopac};
biopac.help    = {'If you choose biopac, physiology records of the subject will be considered as regressors.'};

subj         = cfg_branch;
subj.tag     = 'subj';
subj.name    = 'Subject';
subj.val     = {subj_id age1 anatT1 helmet biopac CWsystem nirs_files protocol TopoData boldmask};%config_path2
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
criugm1.val  = {study_cfg redo1 generic2};
criugm1.prog = @nirs_run_criugm;
criugm1.vout = @nirs_cfg_vout_criugm;
criugm1.help = {'Help'};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_criugm(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
