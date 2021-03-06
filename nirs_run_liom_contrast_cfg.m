function liom_contrast = nirs_run_liom_contrast_cfg
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'Cont');
consess = nirs_spm_get_consess;
display_options = liom_contrast_group_options;

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

StatMethod      = cfg_menu;
StatMethod.tag  = 'StatMethod';
StatMethod.name = 'Statistical method for spatial correlations';
StatMethod.labels = {'EC/LKC','Tube/Bonf'};
StatMethod.values = {1,0};
StatMethod.val  = {1};
StatMethod.help = {'Choose statistical method to account '
    'for false positives due to spatial correlations.'
    '(Family-wise error rate)'
    'Preferred choice: Euler characteristic calculated via Lipschitz-Killing curvature'
    'Other choice: tube formula applicable only to t-statistics for individual sessions,'
    'and Bonferroni correction at patient/subject (group of sessions) level and for F-statistics.'}';

GenerateStats      = cfg_menu;
GenerateStats.tag  = 'GenerateStats';
GenerateStats.name = 'Generate statistics as usual';
GenerateStats.labels = {'Yes','No'};
GenerateStats.values = {1,0};
GenerateStats.val  = {1};
GenerateStats.help = {'This option is for generating interpolated betas only, '
    'it does not generate the usual statistical maps. '
    'These interpolated betas can then be used in the group/anova modules.'}';


UseCorrelRes      = cfg_menu;
UseCorrelRes.tag  = 'UseCorrelRes';
UseCorrelRes.name = 'Use covariance of residuals';
UseCorrelRes.labels = {'Yes','No'};
UseCorrelRes.values = {1,0};
UseCorrelRes.val  = {1};
UseCorrelRes.help = {'This option only applies to the tube formula.'
    'The preferred option is Yes. Select No to obtain the old version'
    'of NIRS_SPM tube calculation, which ignored correlations of residuals'}';

contrast_p_value         = cfg_entry;
contrast_p_value.name    = 'Corrected p_value';
contrast_p_value.tag     = 'contrast_p_value';
contrast_p_value.strtype = 'r';
contrast_p_value.num     = [1 1];
contrast_p_value.val     = {0.05};
contrast_p_value.help    = {'Corrected p_value'};

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


NonlinearEpilepsyOn      = cfg_menu;
NonlinearEpilepsyOn.tag  = 'NonlinearEpilepsyOn';
NonlinearEpilepsyOn.name = 'Automated contrasts for nonlinearities in epilepsy';
NonlinearEpilepsyOn.labels = {'Yes', 'No'};
NonlinearEpilepsyOn.values = {1,0};
NonlinearEpilepsyOn.val    = {0};
NonlinearEpilepsyOn.help = {'This option is ONLY for '
    'Automated contrasts for nonlinearities in epilepsy'}';

TopoData        = cfg_files;
TopoData.tag    = 'TopoData';
TopoData.name   = '(Optional) Select TopoData file';
TopoData.filter = '.mat';
TopoData.num    = [0 1];
TopoData.val{1} = {''};
TopoData.help   = {'This is an option for user to specify a new TopoData,'
    'which is generated during coregistration, instead of the default or current TopoData.'}';

GroupMultiSession           = cfg_menu;
GroupMultiSession.name      = 'Group Multi-Session';
GroupMultiSession.tag       = 'GroupMultiSession';
GroupMultiSession.labels    = {'Yes' 'No'};
GroupMultiSession.values    = {1,0};
GroupMultiSession.val       = {0};
GroupMultiSession.help      = {'Group Multi Session'
    'If selected, with option ProcessContrastBySession set to 0,'
    'Contrasts defined as vectors over all sessions will be treated.'
    'Set to YES, as NO seems to generate bugs (2012-08-09).'}';

user_contrasts         = cfg_branch;
user_contrasts.tag     = 'user_contrasts';
user_contrasts.name    = 'User-defined contrasts'; 
user_contrasts.val     = {GroupMultiSession consess};
user_contrasts.help    = {'User-defined contrasts.'}';

automated_contrasts         = cfg_branch;
automated_contrasts.tag     = 'automated_contrasts';
automated_contrasts.name    = 'Automated contrasts'; 
automated_contrasts.val     = {NonlinearEpilepsyOn};
automated_contrasts.help    = {'Automated contrasts.'}';

ContrastChoice           = cfg_choice;
ContrastChoice.name      = 'Choose contrast choice method';
ContrastChoice.tag       = 'ContrastChoice';
ContrastChoice.values    = {user_contrasts automated_contrasts}; 
ContrastChoice.val       = {user_contrasts}; 
ContrastChoice.help      = {'Choose method to generate contrasts'}'; 

AllowExtrapolation           = cfg_menu;
AllowExtrapolation.name      = 'Allow Extrapolation';
AllowExtrapolation.tag       = 'AllowExtrapolation';
AllowExtrapolation.labels    = {'Yes' 'No'};
AllowExtrapolation.values    = {1,0};
AllowExtrapolation.val       = {1};
AllowExtrapolation.help      = {'Allow Extrapolation.'}';

no_interpolation           = cfg_menu;
no_interpolation.name      = 'Remove interpolation';
no_interpolation.tag       = 'no_interpolation';
no_interpolation.labels    = {'Yes' 'No'};
no_interpolation.values    = {1,0};
no_interpolation.val       = {0};
no_interpolation.help      = {'Remove interpolation entirely.'}';

% ---------------------------------------------------------------------
% sessions Sessions
% ---------------------------------------------------------------------
Sessions         = cfg_entry; %Careful, not the same as sessions in user-defined SPM contrasts
Sessions.name    = 'Sessions';
Sessions.tag     = 'Sessions';
Sessions.strtype = 'r';
Sessions.num     = [0 Inf];
Sessions.val     = {''};
Sessions.help    = {'Enter the numbers of the sessions to include. If not specified, all sessions present in the NIRS matrix will be included.'};

% Executable Branch
liom_contrast      = cfg_exbranch;
liom_contrast.name = 'Liom Contrast Calculations';
liom_contrast.tag  = 'liom_contrast';
%PP removed: liom_contrast_struct
liom_contrast.val  = {NIRSmat redo1 NIRSmatCopyChoice ContrastChoice Sessions ...
    view TopoData GenerateStats StatMethod UseCorrelRes AllowExtrapolation no_interpolation ...
    spatial_LPF contrast_p_value display_options}; %Study_type
liom_contrast.prog = @nirs_run_liom_contrast;
liom_contrast.vout = @nirs_cfg_vout_liom_contrast;
liom_contrast.help = {'Liom Contrast Calculations.'};

function vout = nirs_cfg_vout_liom_contrast(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});