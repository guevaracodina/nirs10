function fcNIRS = nirs_fcNIRS_cfg
% Graphical interface configuration function for functional connectivity mapping
% with near-infrared spectroscopy (fcNIRS) module
%_______________________________________________________________________________
% Copyright (C) 2013 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

fcNIRS          = cfg_choice;
fcNIRS.name     = 'Functional connectivity NIRS (fcNIRS)';
fcNIRS.tag      = 'fcNIRS';
fcNIRS.values   = {nirs_filtdown_cfg nirs_send_email_cfg};
fcNIRS.help     = {'These modules perform resting-state functional connectivity mapping with near-infrared spectroscopy (fcNIRS).'
    'They should be run after the first 3 pre-processing modules, i.e. '
    '1) readBOXY'
    '2) remove_chn_stdev'
    '3) normalize_baseline'
    '4) ODtoHbOHbR'};

% EOF
