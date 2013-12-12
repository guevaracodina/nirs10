function create_roi1 = nirs_create_roi_cfg
% Graphical interface configuration function for nirs_create_roi_run
%_______________________________________________________________________________
% Copyright (C) 2013 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

% Select NIRS.mat
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'ROIs');

% Keep/delete previous ROIs
RemovePreviousROI       = cfg_menu;
RemovePreviousROI.tag   = 'RemovePreviousROI';
RemovePreviousROI.name  = 'Treatment of previous ROIs';
RemovePreviousROI.labels= {'Keep','Remove'};
RemovePreviousROI.values= {false,true};
RemovePreviousROI.val   = {false};
RemovePreviousROI.help  = {'If previous ROIs had been selected, choose whether to keep'
    'or to remove this information.'}';

% Give names for ROIs
select_names            = cfg_menu;
select_names.tag        = 'select_names';
select_names.name       = 'Select names for ROIs';
select_names.labels     = {'Yes','No'};
select_names.values     = {true,false};
select_names.val        = {true};
select_names.help       = {'Option for user to manually enter names of ROIs.'
    'If No is selected, ROI names will be a number (enumeration).'}';
        
ManualROI               = cfg_entry;
ManualROI.tag           = 'ManualROI';
ManualROI.name          = 'Manual ROI selection'; 
ManualROI.strtype       = 'e';
ManualROI.val{1}        = {[1 2]; [3 4 5]};
ManualROI.num           = [Inf 1];
ManualROI.help          = {'Manual ROI selection: specify channels per ROI/seed in a cell column.'}';

ManualROIsel          	= cfg_branch;
ManualROIsel.tag        = 'ManualROIsel';
ManualROIsel.name       = 'Manual ROI selection';
ManualROIsel.val        = {ManualROI};
ManualROIsel.help       = {'Manual ROI selection.'}';


pointNclickROI          = cfg_branch;
pointNclickROI.tag      = 'pointNclickROI';
pointNclickROI.name     = 'Graphic ROI selection: point & click (channels)'; 
pointNclickROI.val      = {};
pointNclickROI.help     = {'Graphic ROI selection: point & click the channels per ROI'
                                }';
                            
AutoROIchoice           = cfg_choice;
AutoROIchoice.name      = 'Choose ROI generation method';
AutoROIchoice.tag       = 'AutoROIchoice';
AutoROIchoice.values    = {ManualROIsel pointNclickROI}; 
AutoROIchoice.val       = {ManualROIsel}; % Default value
AutoROIchoice.help      = {'Choose whether to generate ROI with a list of channels or by clicking'}'; 

% Choose whether interface should ask whether to use a previously saved list of ROI
SelectPreviousROI       = cfg_menu;
SelectPreviousROI.tag   = 'SelectPreviousROI';
SelectPreviousROI.name  = 'Ask for a previous list of ROIs';
SelectPreviousROI.labels= {'No','Yes'};
SelectPreviousROI.values= {0,1};
SelectPreviousROI.val   = {0};
SelectPreviousROI.help  = {'If this option is selected, then before manual selection of ROIs,'
    'The interface will ask the user if they want to use an ROI list from elsewhere.'
    'This ROI list can then be selected by locating the NIRS.mat structure that contains the information.'}';

% Executable Branch
create_roi1             = cfg_exbranch;       % This is the branch that has information about how to run this module
create_roi1.name        = 'Create ROI/seed in channel space';             % The display name
create_roi1.tag         = 'create_roi1'; %Very important: tag is used when calling for execution
create_roi1.val         = {NIRSmat redo1 NIRSmatCopyChoice RemovePreviousROI  ...
    select_names AutoROIchoice SelectPreviousROI };    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
create_roi1.prog        = @nirs_create_roi_run;  % A function handle that will be called with the harvested job to run the computation
create_roi1.vout        = @nirs_cfg_vout_create_roi; % A function handle that will be called with the harvested job to determine virtual outputs
create_roi1.help        = {'Create regions of interest/seeds.'};
return

%make PAT.mat available as a dependency
function vout = nirs_cfg_vout_create_roi(job)
vout                    = cfg_dep;                  % The dependency object
vout.sname              = 'NIRS.mat';               % Displayed dependency name
vout.src_output         = substruct('.','NIRSmat');
vout.tgt_spec           = cfg_findspec({{'filter','mat','strtype','e'}});

% EOF
