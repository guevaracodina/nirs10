function extract_map_data = nirs_run_extract_map_data_cfg

NIRSmat         = cfg_files; %Select NIRS.mat for this subject
NIRSmat.name    = 'NIRS.mat'; % The displayed name
NIRSmat.tag     = 'NIRSmat';       %file names
NIRSmat.filter  = 'mat';
NIRSmat.ufilter = '^NIRS.mat$';
NIRSmat.num     = [1 Inf];     % Number of inputs required
NIRSmat.help    = {'Select NIRS.mat for the subject(s).'}; % help text displayed



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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Extract map data from group and session analyses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

extract_TOPO_map         = cfg_files;
extract_TOPO_map.name    = 'Select TOPO.mat for map to use';
extract_TOPO_map.tag     = 'extract_TOPO_map';
extract_TOPO_map.filter = 'mat';
extract_TOPO_map.ufilter = '^TOPO.mat$';
extract_TOPO_map.num     = [1 1];
extract_TOPO_map.val     = {''};
extract_TOPO_map.help    = {'Optional. Default is TOPO.mat saved in specified NIRS.mat.'
    'This is the TOPO.mat that will be used to select pixels of interest to extract'}';

extract_TOPO_con         = cfg_files;
extract_TOPO_con.name    = 'Select (big_)TOPO.mat for contrasts to use';
extract_TOPO_con.tag     = 'extract_TOPO_con';
extract_TOPO_con.filter = 'mat';
extract_TOPO_con.ufilter = 'TOPO.mat$';
extract_TOPO_con.num     = [1 1];
extract_TOPO_con.val     = {''};
extract_TOPO_con.help    = {'Optional. Default is TOPO.mat saved in specified NIRS.mat.'
    'This is the TOPO.mat that will be used to get the betas associated with specific contrasts.'
    'For a multi-subject study, this would the big_TOPO.mat file'}';

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

% Executable Branch
extract_map_data      = cfg_exbranch;
extract_map_data.name = 'Extract map data';
extract_map_data.tag  = 'extract_map_data';
extract_map_data.val  = {NIRSmat extract_TOPO_map extract_TOPO_con view extract_base_contrast ...
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