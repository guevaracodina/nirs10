function bases = nirs_dfg_hrf
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
% ---------------------------------------------------------------------
% hrf Canonical HRF
% ---------------------------------------------------------------------
hrf         = cfg_branch;
hrf.tag     = 'hrf';
hrf.name    = 'Canonical HRF';
hrf.val     = {derivs };
hrf.help    = {'Canonical Hemodynamic Response Function. This is the default option. Contrasts of these effects have a physical interpretation and represent a parsimonious way of characterising event-related responses. This option is also useful if you wish to look separately at activations and deactivations (this is implemented using a t-contrast with a +1 or -1 entry over the canonical regressor). '};
% ---------------------------------------------------------------------
% eHRF: epilepsy HRF -- precise form to be studied and determined --
% for now: this is an HRF that falls off very slowly
% ---------------------------------------------------------------------

eHRF         = cfg_branch;
eHRF.tag     = 'eHRF';
eHRF.name    = 'Epilepsy HRF';
eHRF.val     = {derivs };
eHRF.help    = {'Hemodynamic Response Function for epilepsy. '
    'For now: this is an HRF that falls off very slowly and has no undershoot'}';

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
bases.values  = {hrf gamma eHRF}; %{hrf fourier fourier_han gamma fir };
