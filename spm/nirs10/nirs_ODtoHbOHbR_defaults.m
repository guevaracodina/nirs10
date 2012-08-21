function nirs_ODtoHbOHbR_defaults

global nirs10

%nirs10.preprocessNIRS.ODtoHbOHbR.Normalize_OD = 0; %0: Median; 1: Initial value 2: Mean
%nirs10.ODtoHbOHbR.subject_age = 25; %Age in years
nirs10.preprocessNIRS.ODtoHbOHbR.PVF = [50 50];
nirs10.preprocessNIRS.ODtoHbOHbR.DPF = [6.51 5.86]; % FOR 690 and 832 nm
% From Duncan et al., Phys. Med. Biol. 1995 40(2):295-304.
%nirs10.preprocessNIRS.ODtoHbOHbR.threshold = 0.1;
end