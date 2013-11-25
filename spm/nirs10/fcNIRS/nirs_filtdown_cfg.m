function filtdown1 = nirs_filtdown_cfg
% Graphical interface configuration function for ioi_filtdown_run
%_______________________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

% Select IOI.mat
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'filtDown');
% % Choose ROI selection method (all/selected)
% ROI_choice          = ioi_dfg_ROI_choice;
% % Choose session selection method (all/selected)
% session_choice      = ioi_dfg_session_choice;
% % Colors to include (OD,HbO,HbR,HbT,Flow)
% IC                  = ioi_dfg_include_colors(0,1,1,1,1);
% % Bandpass filtering
% bpf                 = ioi_bpf_cfg(1, [0.009 0.08], 4, 'butter');
% % Generate / save figures
% [generate_figures ...
%     save_figures]       = ioi_dfg_generate_figures;

% Downsampling frequency
downFreq            = cfg_entry;
downFreq.name       = 'Downsampling frequency in Hz';   % The displayed name
downFreq.tag        = 'downFreq';                       % file names
downFreq.strtype    = 'r';                              % Real numbers
downFreq.num        = [1 1];                            % Number of inputs required
downFreq.val        = {1};                              % Default value
downFreq.help       = {'Enter downsampling frequency in Hz. A target sampling frequency will be generated, which may however be only approximately equal to the specified downsampling frequency, but it will correspond to the actual frequency of selecting every Nth point'};

% Filter and downsample whole images
wholeImage          = cfg_menu;
wholeImage.tag      = 'wholeImage';
wholeImage.name     = 'Filter and downsample whole image time-series';
wholeImage.labels   = {'Yes','No'};
wholeImage.values   = {1,0};
wholeImage.val      = {1};
wholeImage.help     = {'Filter and downsample whole image time-series. It creates a new sub-folder for each session'};

% Executable Branch
filtdown1           = cfg_exbranch; % This is the branch that has information about how to run this module
filtdown1.name      = 'Temporal filtering & downsampling';             % The display name
filtdown1.tag       = 'filtdown1'; %Very important: tag is used when calling for execution
filtdown1.val       = {NIRSmat redo1 NIRSmatCopyChoice downFreq wholeImage};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
filtdown1.prog      = @nirs_filtdown_run;  % A function handle that will be called with the harvested job to run the computation
filtdown1.vout      = @nirs_cfg_vout_filtdown; % A function handle that will be called with the harvested job to determine virtual outputs
filtdown1.help      = {'Temporal band-pass filtering and downsampling of a given time trace [HbO/HbR], either on a seed or on all the channels.'};

return

% Make IOI.mat available as a dependency
function vout = nirs_cfg_vout_filtdown(job)
vout                = cfg_dep; % The dependency object
vout.sname          = 'NIRS.mat'; % Displayed dependency name
vout.src_output     = substruct('.','NIRSmat'); %{1}); %,'IOImat');
vout.tgt_spec       = cfg_findspec({{'filter','mat','strtype','e'}});

% EOF

