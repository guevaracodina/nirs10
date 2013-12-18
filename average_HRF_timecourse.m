function out = average_HRF_timecourse(Y, onsets, fs, time_before_spk, time_after_spk, EvtofInterest, chromophore, avoidance)

Nc = size(Y,1);
switch chromophore
    case 0 %HbO
        YC = Y(1:(Nc/3),:);
    case 1 %HbR
        YC = Y((Nc/3+1):(2*Nc/3),:);
    case 2 %HbT
        YC = Y((2*Nc/3+1):Nc,:);
end

spk_ons = onsets(EvtofInterest).ons;
disp(['Total number of ' int2str(length(spk_ons)) ' Event of interest in the session']);

% Perform avoidance
if avoidance.enabled
    EvtofNonInterest = avoidance.EvtofNonInterest;
    for e0 = 1 : length(EvtofNonInterest)
        event = EvtofNonInterest(e0);
        EvtofNonInterest_ons = onsets(event).ons;
        EvtofNonInterest_dur = onsets(event).dur;
        EvtofNonInterest_offs = EvtofNonInterest_ons + EvtofNonInterest_dur;
        spk_ons_filtered = [];
        for i0 = 1 : length(spk_ons)
            spk_timing = spk_ons(i0);
            EvtofNonInterest_late = find(EvtofNonInterest_ons >= spk_timing);
            if ~isempty(EvtofNonInterest_late)
                EvtofNonInterest_nextonset= EvtofNonInterest_late(1);
                EvtofNonInterest_nextonsettime = EvtofNonInterest_ons(EvtofNonInterest_nextonset);
            else
                EvtofNonInterest_nextonsettime = Inf;
            end
            EvtofNonInterest_previous = find(EvtofNonInterest_offs <= spk_timing);
            if ~isempty(EvtofNonInterest_previous)
                EvtofNonInterest_previousonset = EvtofNonInterest_previous(end);
                EvtofNonInterest_previousonsettime = EvtofNonInterest_offs(EvtofNonInterest_previousonset);
            else
                EvtofNonInterest_previousonsettime = -Inf;
            end
            
            if (spk_timing >= EvtofNonInterest_previousonsettime - time_before_spk) && (spk_timing <= EvtofNonInterest_nextonsettime - time_after_spk)
                spk_ons_filtered = [spk_ons_filtered spk_timing];
            end
        end
        spk_ons = spk_ons_filtered;
    end
    disp(['After avoidance, total number of ' int2str(length(spk_ons)) ' Event of interest left']);
end

%Discard events that are too close to the start or to the end
spk_NotTooEarly = find(spk_ons >= -time_before_spk);
if ~isempty(spk_NotTooEarly)
    spk_ons = spk_ons(spk_NotTooEarly);
end
spk_NotTooLate = find(spk_ons <= (size(YC,2)/fs - time_after_spk));
if ~isempty(spk_NotTooLate)
    spk_ons = spk_ons(spk_NotTooLate);
end
RestSpkNum = length(spk_ons);
disp([' Final Spk Num for averaging:' int2str(RestSpkNum)]);

%Average
PointsBefore = round(abs(time_before_spk)*fs);
PointsAfter = round(abs(time_after_spk)*fs);
total_num = PointsBefore+PointsAfter+1;
HRF_all = zeros(Nc/3,total_num);

for i0 = 1 : RestSpkNum
    spk_timing = round(spk_ons(i0)*fs);
    HRF_all = HRF_all + YC(:,(spk_timing - PointsBefore) : (spk_timing + PointsAfter));
end
HRF_avg = HRF_all/RestSpkNum;
out = HRF_avg;

