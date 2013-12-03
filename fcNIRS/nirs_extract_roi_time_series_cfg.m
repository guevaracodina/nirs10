function extract_roi1 = nirs_extract_roi_time_series_cfg
% Graphical interface configuration function for nirs_extract_roi_time_series_run
%_______________________________________________________________________________
% Copyright (C) 2013 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

% Select NIRS.mat
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'ROIs');
% Choose ROI selection method (all/selected)
% ROI_choice              = pat_roi_choice_cfg;
% Colors to include (HbT, SO2, Bmode)
% IC                      = pat_include_colors_cfg(true, true);


% Extract mean signal (from the brain mask pixels)
extractBrainMask        = cfg_menu;
extractBrainMask.tag    = 'extractBrainMask';
extractBrainMask.name   = 'Extract brain mask signal';
extractBrainMask.labels = {'Yes','No'};
extractBrainMask.values = {true, false};
extractBrainMask.val    = {false};   % Default value = 0
extractBrainMask.help   = {'Extract mean signal from the non-masked brain pixels'};
% Extract mean signal (from the brain mask pixels) -- this time using an
% activation map
% activMask_choice        = pat_activation_mask_choice_cfg(false);

% Executable Branch
extract_roi1            = cfg_exbranch;       % This is the branch that has information about how to run this module
extract_roi1.name       = 'Extract ROI/seed in channel space';             % The display name
extract_roi1.tag        = 'extract_roi1'; %Very important: tag is used when calling for execution
extract_roi1.val        = {NIRSmat redo1 NIRSmatCopyChoice ...
    extractBrainMask  };    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
extract_roi1.prog       = @nirs_extract_roi_time_series_run;  % A function handle that will be called with the harvested job to run the computation
extract_roi1.vout       = @nirs_cfg_vout_extract_roi; % A function handle that will be called with the harvested job to determine virtual outputs
extract_roi1.help       = {'Create regions of interest.'};

return

%make PAT.mat available as a dependency
function vout = nirs_cfg_vout_extract_roi(job)
vout                    = cfg_dep;                  % The dependency object
vout.sname              = 'NIRS.mat';               % Displayed dependency name
vout.src_output         = substruct('.','NIRSmat');
vout.tgt_spec           = cfg_findspec({{'filter','mat','strtype','e'}});

% EOF
