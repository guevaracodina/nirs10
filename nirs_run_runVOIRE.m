function out = nirs_run_runVOIRE(job)
%__________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire
% Cl�ment Bonn�ry
% 2010-12

% le rythme cardiaque me donne
%Friston_nonlinear_2000 : rCBF lineaire avec sunaptic activity
% CBV grossi�rement �gal � HbT : diff�rence est dans la plasme qui est
% l'autre partie du sang !!
% Cox : oxyg�nation du sang : grossi�rement HbO/HbT (on tient pas compte du plasma...)
% D�bit cardiaque = freq card * volume �ject�
prefix = 'h';

if job.heart_pace==1
    jobHP.NIRSmat = job.NIRSmat;
    jobHP.STFT_param.win_type = 0;
    jobHP.STFT_param.win_width =6;
    jobHP.STFT_param.Nprobe = 200;
    jobHP.STFT_param.fft_size = 1024;
    
    outHP = nirs_run_criugm_paces(jobHP);
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
%         [dir1 fil1 ext1] = fileparts(rDtp{f});
%         if strcmp(ext1,'.nirs')
%             try load(fullfile(dir1,[prefix,fil1,ext1]),'-mat'); 
%                 c=c';
%             catch
%                 disp('chargement �chou�');
%             end
%         end
        
        Cid = NIRS.Cf.H.C.id;
        COx = zeros(size(c,1),size(c,2)/2);
        for iC = 1:size(Cid,2)/2
             COx(:,iC) = c(iC)/(c(iC)+c(iC+size(c,2)/2));
        end
        save(fullfile(NIRS.Dt.s.p,['Cox' fil1 '.mat']),'COx');
    end
    
catch
    disp('Could not evaluate COx');
end
out.NIRSmat = job.NIRSmat;