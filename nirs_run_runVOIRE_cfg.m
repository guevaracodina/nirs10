function runVOIRE1 = nirs_run_runVOIRE_cfg
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'voire');

%%%%%%%%%%%%%%%

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
criugm_paces1.val  = {NIRSmat NIRSmatCopyChoice heart_rate_cfg ...
    remove_no_heartbeat save_heart_rate_figure display_heart_rate_figure};
criugm_paces1.prog = @nirs_run_criugm_paces;
criugm_paces1.vout = @nirs_cfg_vout_criugm_paces;
criugm_paces1.help = {['Preprocessing step: Extract heart rate and, if desired, ',...
    'remove channels without a clear detectable heart rate.']}';

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
runVOIRE1.val  = {NIRSmat redo1 NIRSmatCopyChoice heart_pace criugm_paces1};
runVOIRE1.prog = @nirs_run_runVOIRE;
runVOIRE1.vout = @nirs_cfg_vout_runVOIRE;
runVOIRE1.help = {'.'};

function vout = nirs_cfg_vout_runVOIRE(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_criugm_paces(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});