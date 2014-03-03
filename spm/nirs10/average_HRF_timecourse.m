function [HRF_avg HRF_std] = average_HRF_timecourse(Y, onsets, fs, time_before_spk, time_after_spk, EvtofInterest, EvtofInterest_name, chromophore, avoidance, sub)

Nc = size(Y,1);
switch chromophore
    case 1 %HbO
        YC = Y(1:(Nc/3),:);
    case 2 %HbR
        YC = Y((Nc/3+1):(2*Nc/3),:);
    case 3 %HbT
        YC = Y((2*Nc/3+1):Nc,:);
end

ons_sz = [];
ons_spk = [];
for i0 = 1 : length(onsets)
    name_ons = onsets(i0).name{1};
    if ~isempty(EvtofInterest_name) && strcmpi(name_ons, EvtofInterest_name)
        spk_ons = onsets(i0).ons;
        EvtofInterest_name = name_ons;
    elseif  strcmpi(name_ons, 'sz')
        ons_sz = [ons_sz i0];
    elseif (~isempty(strfind(name_ons, 'spk'))) || (~isempty(strfind(name_ons, 'Spk'))) || ~isempty(strfind(name_ons, 'SPK'))
        if i0 ~= EvtofInterest
            ons_spk = [ons_spk i0];
        end
    end
end
if EvtofInterest ~= 0
    spk_ons = onsets(EvtofInterest).ons;
    EvtofInterest_name = onsets(EvtofInterest).name;
end
disp(['Total number of ' int2str(length(spk_ons)) ' Event of interest (' EvtofInterest_name ') in the session']);
disp('-----------------------------------------------------------------');

%Discard events that are too close to the start or to the end
spk_NotTooEarly = find(spk_ons >= -time_before_spk);
if ~isempty(spk_NotTooEarly)
    spk_ons = spk_ons(spk_NotTooEarly);
end
spk_NotTooLate = find(spk_ons <= (size(YC,2)/fs - time_after_spk));
if ~isempty(spk_NotTooLate)
    spk_ons = spk_ons(spk_NotTooLate);
end

% Perform avoidance
if avoidance.enabled
    if avoidance.type == 1
        EvtofNonInterest_sz = avoidance.EvtofNonInterest_sz;
        EvtofNonInterest_spk = avoidance.EvtofNonInterest_spk;
        if (length(EvtofNonInterest_sz) == 1) && (EvtofNonInterest_sz == 0)
            EvtofNonInterest_sz = ons_sz;
        end
        if (length(EvtofNonInterest_spk) == 1) && (EvtofNonInterest_spk == 0)
            EvtofNonInterest_spk = ons_spk;
        end
        if length(onsets) > 1
            show_ons_av(onsets, EvtofNonInterest_sz, EvtofNonInterest_spk);
            spk_ons = EvtofNonInterest_avoid(onsets, spk_ons, EvtofNonInterest_sz, EvtofNonInterest_spk, time_before_spk, time_after_spk);
        end
        spk_valid = spk_ons;
    end
    if avoidance.type == 2
        spk_valid = spk_ons;
        spk_ons = Internal_avoid(spk_ons, time_before_spk, time_after_spk);
    end
    if avoidance.type == 3
        EvtofNonInterest_sz = avoidance.EvtofNonInterest_sz;
        EvtofNonInterest_spk = avoidance.EvtofNonInterest_spk;
        if (length(EvtofNonInterest_sz) == 1) && (EvtofNonInterest_sz == 0)
            EvtofNonInterest_sz = ons_sz;
        end
        if (length(EvtofNonInterest_spk) == 1) && (EvtofNonInterest_spk == 0)
            EvtofNonInterest_spk = ons_spk;
        end        
        if length(onsets) > 1
            show_ons_av(onsets, EvtofNonInterest_sz, EvtofNonInterest_spk);
            spk_ons = EvtofNonInterest_avoid(onsets, spk_ons, EvtofNonInterest_sz, EvtofNonInterest_spk, time_before_spk, time_after_spk);
        end
        disp(['After avoidance of event of non-interest, total number of ' int2str(length(spk_ons)) ' Event of interest left']);
        spk_valid = spk_ons;
        spk_ons = Internal_avoid(spk_ons, time_before_spk, time_after_spk);
        disp(['After avoidance of interaction time-courses (internal), total number of ' int2str(length(spk_ons)) ' Event of interest left']);
    end
    disp(['After avoidance, total number of ' int2str(length(spk_ons)) ' Event of interest left']);
end

RestSpkNum = length(spk_ons);
disp('-----------------------------------------------------------------');
disp(['Final Spk Num for averaging:' int2str(RestSpkNum)]);

%Average
PointsBefore = round(abs(time_before_spk)*fs);
PointsAfter = round(abs(time_after_spk)*fs);
total_num = PointsBefore+PointsAfter+1;
HRF_all = zeros(Nc/3,total_num);
HRF_iv = zeros(total_num * (Nc/3), RestSpkNum);

for i0 = 1 : RestSpkNum
    spk_timing = round(spk_ons(i0)*fs);
    HRF_id = YC(:,(spk_timing - PointsBefore) : (spk_timing + PointsAfter));
    HRF_iv(:,i0) = reshape(HRF_id', total_num * (Nc/3), 1);
    HRF_all = HRF_all + YC(:,(spk_timing - PointsBefore) : (spk_timing + PointsAfter));
end
HRF_avg = HRF_all/RestSpkNum;
%HRF_iv = HRF_iv';
%HRF_std = std(HRF_iv,0,2); %standard deviation
HRF_std = std(HRF_iv,0,2)/sqrt(RestSpkNum); %standard error

if ~sub %No subtraction
else
    
end

function spk_ons = EvtofNonInterest_avoid(onsets, spk_ons, EvtofNonInterest_sz, EvtofNonInterest_spk, time_before_spk, time_after_spk)
%Sz-like events
EvtofNonInterest_ons_sz = [];
EvtofNonInterest_offs_sz = [];
for e0 = 1 : length(EvtofNonInterest_sz)
    event = EvtofNonInterest_sz(e0);
    if event ~= 0
        EvtofNonInterest_ons0 = onsets(event).ons;
        EvtofNonInterest_dur0 = onsets(event).dur;
        EvtofNonInterest_offs0 = EvtofNonInterest_ons0 + EvtofNonInterest_dur0;
        EvtofNonInterest_ons_sz = [EvtofNonInterest_ons_sz EvtofNonInterest_ons0];
        EvtofNonInterest_offs_sz = [EvtofNonInterest_offs_sz EvtofNonInterest_offs0];
    end
end
%Spk-like events
EvtofNonInterest_ons_spk = [];
EvtofNonInterest_offs_spk = [];
for e0 = 1 : length(EvtofNonInterest_spk)
    event = EvtofNonInterest_spk(e0);
    if event ~= 0
        EvtofNonInterest_ons0 = onsets(event).ons;
        EvtofNonInterest_offs0 = EvtofNonInterest_ons0; %zero-duration
        EvtofNonInterest_ons_spk = [EvtofNonInterest_ons_spk EvtofNonInterest_ons0];
        EvtofNonInterest_offs_spk = [EvtofNonInterest_offs_spk EvtofNonInterest_offs0];
    end
end
%Combine
EvtofNonInterest_ons = [EvtofNonInterest_ons_sz EvtofNonInterest_ons_spk];
EvtofNonInterest_offs = [EvtofNonInterest_offs_sz EvtofNonInterest_offs_spk];
%Perform filtering
spk_ons = EvtofNonInterest_filter(spk_ons, EvtofNonInterest_ons, EvtofNonInterest_offs, time_before_spk, time_after_spk);

function spk_ons_filtered = EvtofNonInterest_filter(spk_ons, EvtofNonInterest_ons, EvtofNonInterest_offs, time_before_spk, time_after_spk)
spk_tmp = spk_ons;
for i0 = 1 : length(EvtofNonInterest_ons)
    time_before = EvtofNonInterest_ons(i0) - time_after_spk + time_before_spk;
    time_after = EvtofNonInterest_offs(i0) - time_before_spk + time_after_spk;
    spk_idx_tmp1 = find(spk_tmp >= time_after);
    spk_idx_tmp2 = find(spk_tmp <= time_before);
    spk_idx = [spk_idx_tmp2 spk_idx_tmp1];
    spk_tmp = spk_tmp(spk_idx);
    if isempty(spk_tmp)
        break;
    end
end
spk_ons_filtered = spk_tmp;

function spk_ons = Internal_avoid(spk_ons, time_before_spk, time_after_spk)
if length(spk_ons) > 2
    flag = ones(1,length(spk_ons));
    spk_ons_end = spk_ons + time_after_spk;
    spk_ons_start = spk_ons + time_before_spk;
    spk_ons_end_shifted = circshift(spk_ons_end,[0,1]);
    spk_ons_start_shifted = circshift(spk_ons_start,[0,-1]);
    spk_timing_diff1 = spk_ons_start - spk_ons_end_shifted;
    spk_timing_diff2 = spk_ons_end - spk_ons_start_shifted;
    spk_idx1 = find(spk_timing_diff1 < 0);
    flag(spk_idx1) = 0;
    spk_idx2 = find(spk_timing_diff2 > 0);
    flag(spk_idx2) = 0;
    if spk_timing_diff1(2) >= 0
        flag(1) = 1;
    end
    if spk_timing_diff2(end-1) <= 0
        flag(end) = 1;
    end
    spk_ons = spk_ons(find(flag));
elseif length(spk_ons) == 2
    if (spk_ons(2) - spk_ons(1)) < (time_after_spk - time_before_spk)
        spk_ons = [];
    end
end

function show_ons_av(ons, sz, spk)
disp('Onsets avoided:');
for i0 = 1 : length(sz)
    disp('Sz-like events:');
    idx = sz(i0);
    name = ons(idx).name{1};
    n = length(ons(idx).ons);
    disp([name ' ' int2str(n)]);
end
for i0 = 1 : length(spk)
    disp('Spk-like events:');
    idx = spk(i0);
    name = ons(idx).name{1};
    n = length(ons(idx).ons);
    disp([name ' ' int2str(n)]);
end
    
    



    %         spk_ons_filtered = [];
    %         for i0 = 1 : length(spk_ons)
    %             spk_timing = spk_ons(i0);
    %             EvtofNonInterest_late = find(EvtofNonInterest_ons >= spk_timing);
    %             if ~isempty(EvtofNonInterest_late)
    %                 EvtofNonInterest_nextonset= EvtofNonInterest_late(1);
    %                 EvtofNonInterest_nextonsettime = EvtofNonInterest_ons(EvtofNonInterest_nextonset);
    %             else
    %                 EvtofNonInterest_nextonsettime = Inf;
    %             end
    %             EvtofNonInterest_previous = find(EvtofNonInterest_offs <= spk_timing);
    %             if ~isempty(EvtofNonInterest_previous)
    %                 EvtofNonInterest_previousonset = EvtofNonInterest_previous(end);
    %                 EvtofNonInterest_previousonsettime = EvtofNonInterest_offs(EvtofNonInterest_previousonset);
    %             else
    %                 EvtofNonInterest_previousonsettime = -Inf;
    %             end
    %             
    %             if (spk_timing >= EvtofNonInterest_previousonsettime - time_before_spk) && (spk_timing <= EvtofNonInterest_nextonsettime - time_after_spk)
    %                 spk_ons_filtered = [spk_ons_filtered spk_timing];
    %             end
    %         end
    %         spk_ons = spk_ons_filtered;