function d = nirs_new_remove_jumps(d)
%[ dData, minmax, stats ] = AnalyzeEdges( data, scales, thresholds, timestep, startTime, endTime, tranRad );
try
    N = size(d,2);
    wi = 5;
    level_after = 0;
    level_before = 0;
    for i0 = 1:size(d,1)
        %test
        %d(i0,1:1000) = d(i0,1:1000)-1e4;
        td = d(i0,:);
        [dData{i0},minmax{i0},stats{i0}] = nirs_AnalyzeEdges(td',[1],[0.99]);
        mm= find(minmax{i0});
        for j0 = 1:(length(mm)+1)
            if j0 == 1
                si = 1;
            else
                si = min(mm(j0-1)+1,N);
            end
            if j0 == length(mm)+1
                ei = N;
            else
                ei = mm(j0);
            end
            
            Mwi = max(1,ei-wi); %catch case where ei is small
            td(si:Mwi) = td(si:Mwi) - level_after+level_before;
            
            Bwi = min(N,ei+wi);
            B2wi = min(N,ei+2*wi);
            M2wi = max(1,ei-2*wi);
            if j0 < length(mm)+1
                level_after = mean(td(Bwi:B2wi));
                level_before = mean(td(M2wi:Mwi));
            end
            td(Mwi:ei) = level_before;
        end
        %figure; plot(td); hold on; plot(d(i0,:),'r')
        d(i0,:) = td;
    end
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
end