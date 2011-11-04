function out = offset_nirs(job)
load(job.NIRSmat);
iR = job.iR;
offset = 0;

% loop over sessions
% for iR=1:size(job.subj(1,is).nirs_files,1)
    subj.nirs = load(NIRS.Dt.fir.pp(1,1).p{iR,1}, '-mat'); % get the IO data
    aux = subj.nirs.aux;
    t = subj.nirs.t;
    % Retrieve 4 channels.
    aux_n = aux(:,1:2);  % conditions are coded over 2 first bits of auxilaries
    
    % Normalize auxiliaries.
    aux_n = aux_n - min(min(aux_n, [], 1));
    aux_n = aux_n / max(max(aux_n, [], 1));
    seuil = 0.8*max(aux_n);
    
    % Attribute 0 and 1 to arbitrary values of aux.nirs
    trig = zeros(size(t,1),size(aux_n,2));
    Aux1 = [];
    Aux2 = [];
    
    for it=1:size(t)
        for iaux=1:size(aux_n,2)
            if aux_n(it,iaux) > seuil(iaux)
                trig(it,iaux) = 1;
            else trig(it,iaux) = 0;
            end
        end
    end
    
    Aux1 = trig(:,1);
    Aux2 = trig(:,2);
    
    Offset1 = find(Aux1,1,'first');
    Offset2 = find(Aux2,1,'first');
    offset = min(Offset1, Offset2)*0.04;  
    
%     clear subj.nirs aux_n
    
% end

out = offset;