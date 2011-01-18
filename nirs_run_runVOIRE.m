function out = nirs_run_runVOIRE(job)
%__________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire
% Clément Bonnéry
% 2010-12

% le rythme cardiaque me donne
%Friston_nonlinear_2000 : rCBF lineaire avec sunaptic activity
% CBV grossièrement égal à HbT : différence est dans la plasme qui est
% l'autre partie du sang !!
% Cox : oxygénation du sang : grossièrement HbO/HbT (on tient pas compte du plasma...)
% Débit cardiaque = freq card * volume éjecté

% changes in the concentration of total Hb as a measure for changes in cerebral blood volume (CBV)
% and hemoglobin difference (HbD; oxy- minus deoxy-Hb) : CBF (Tsuji, M., Duplessis, A., Taylor, G., Crocker, R., Volpe, J.J., 1998. Near infrared spectroscopy detects cerebral ischemia during hypotension in piglets. Pediatr. Res. 44, 591–595.)

prefix = 'h';

if job.heart_pace==1
%     jobHP.NIRSmat = job.NIRSmat;
%     jobHP.STFT_param.win_type = 0;
%     jobHP.STFT_param.win_width =6;
%     jobHP.STFT_param.Nprobe = 200;
%     jobHP.STFT_param.fft_size = 1024;
%     remove_no_heartbeat 
%     detect_wavelength 
%     MinHeartRate 
%     MaxHeartRate
%     InternalMinHeartRate 
%     InternalMaxHeartRate 
%     MaxHeartStdev
    
    outHP = nirs_run_paces(job.paces1);
    hp = outHP.heartpace;
end

load(job.NIRSmat{:});

jobHb.age = NIRS.Dt.s.age;
jobHb.NIRSmat = job.NIRSmat;
jobHb.Normalize_OD = 0;
jobHb.subject_age = NIRS.Dt.s.age;
jobHb.PVF = [50;50];
jobHb.threshold = 0.1;
nirs_lpf2.lpf_gauss2.fwhm1 =1.5;
nirs_lpf2.lpf_gauss2.downsamplingFactor =10;
nirs_lpf2.lpf_gauss2.downsampleWhen = 1;
jobHb.nirs_lpf2 = nirs_lpf2;
outHb = nirs_run_ODtoHbOHbR(jobHb);

%use last step of preprocessing
lst = length(NIRS.Dt.fir.pp);
rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
try
    %loop over data files
    for f=1:size(rDtp,1)
        c = fopen_NIR(rDtp{f});
        [~,fil1,~] = fileparts(rDtp{f});
        
        Cid = NIRS.Cf.H.C.id;
        COx = zeros(size(c,1)/2,size(c,2)); % Saturation O2
        CBF = zeros(size(c,1)/2,size(c,2)); % CBF
        HbT = zeros(size(c,1)/2,size(c,2)); % HbT/CBV
        
        for iC = 1:size(Cid,2)/2
            CBF(iC,:) = (c(iC,:)-c(iC+size(c,1)/2,:)); % HbO - HbR : CBF
            HbT(iC,:) = (c(iC,:)+c(iC+size(c,1)/2,:)); % Hb total : CBV
            COx(iC,:) = c(iC,:)./HbT(iC,:);
        end
        save(fullfile(NIRS.Dt.s.p,['Cox' fil1 '.mat']),'COx');
        save(fullfile(NIRS.Dt.s.p,['CBF' fil1 '.mat']),'CBF');
        save(fullfile(NIRS.Dt.s.p,['HbT' fil1 '.mat']),'HbT');
        % on ajoute aussi à la matrice NIRS
    end
    
catch
    disp('Could not evaluate COx or HbT or CBF');
end
out.NIRSmat = job.NIRSmat;