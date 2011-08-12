function out = nirs_run_runOvertraining(job)
% fin analyse Olivier -- Overtraining
try
    viewer_ON = job.View6Projections;
catch
    viewer_ON = 0;
end
try
    NewNIRSdir = job.NewDirCopyNIRS.CreateNIRSCopy.NewNIRSdir;
    NewDirCopyNIRS = 1;
catch
    NewDirCopyNIRS = 0;
end

% Loop over subjects
for iSubj=1:size(job.NIRSmat,1)
    
    % Load NIRS.mat
    try
        NIRS = [];
        load(job.NIRSmat{iSubj,1});
        
%         lst = size(NIRS.Dt.fir.pp,2);
%         p = NIRS.Dt.fir.pp(1,lst).p;
        
        for iSess =1:size(NIRS.Dt.fir.Sess,2)
            % on moyenne aux 10sec
%             D = load(p{iSess,1},'-mat');
            D_time = load(NIRS.Dt.fir.pp(1,1).p{iSess,1},'-mat');
            [dir,n,e]= fileparts(NIRS.Dt.fir.pp(1,1).p{iSess,1});
            
            ind_end = (max(D_time.t)-20)*25+1;
            ind_start =20*25+1;
            
            names        = {'baseline',['cond_' n(end)]};
            
            onsets{1}    = [1,ind_end];
            durations{1} = [20,20];
            
            onsets{2}    = ind_start;
            durations{2} = 540;
            
            save(fullfile(dir,['onsets_' n]),'onsets','names','durations');
        end
        
        if NewDirCopyNIRS
            [dirN fil1 ext1] =fileparts(job.NIRSmat{iSubj,1});
            dir2 = [dirN filesep NewNIRSdir];
            if ~exist(dir2,'dir'), mkdir(dir2); end;
            newNIRSlocation = fullfile(dir2,'NIRS.mat');
            save(newNIRSlocation,'NIRS');
            job.NIRSmat{iSubj,1} = newNIRSlocation;
        else
            save(job.NIRSmat{iSubj,1},'NIRS');
        end
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        disp(['Analysis Overtraining failed for the ' int2str(iSubj) 'th subject.']);
    end
end
out.NIRSmat = job.NIRSmat;%job.NIRSmat{iSubj};


%         for iSess =1:size(NIRS.Dt.fir.Sess,2)
%             % on moyenne aux 10sec
%             D = load(p{iSess,1},'-mat');
%             D_time = load(NIRS.Dt.fir.pp(1,1).p{iSess,1},'-mat');
%             [dir,n,e]= fileparts(NIRS.Dt.fir.pp(1,1).p{iSess,1});
%             
%                             count =0;
%                 d_moyenne =zeros(251,56);
%                 
%             try
%                 %20sec de baseline
%                 ind_start = 20*25+1;
%                 %tache de 600-40 = 560sec soit 9min et 20sec
%                 %20sec de baseline
%                 ind_end = 580*25+1;
%                 test = D.d(ind_end,1);
%             catch
%                 ind_end = (max(D_time.t)-20)*25+1;
%                 ind_start = ind_end - 540*25+1;
%             end
%             for iSlap=ind_start:10*25+1:ind_end-(10*25+1)
%                 d_moyenne = d_moyenne+D.d(iSlap:iSlap+10*25,:);
%                 count = count+1;
%             end
%             d_moyenne = d_moyenne/count;
%             save(fullfile(dir,['M10s_' n '.mat']),'d_moyenne');
%         end