function [bpi bpd] = find_marker_durations(bpi)
gbp = [];
bpd = [];
if bpi
    count = 1;
    if length(bpi)>1
        for i=1:length(bpi)-1        
            if bpi(i+1)-bpi(i) > 1
                %there is now a gap
                gbp = [gbp bpi(i)-count+1]; %take left-most point
                bpd = [bpd count];
                count = 1;
            else
                count = count + 1;
            end    
        end
        %add last marker
        if bpi(end)-bpi(end-1) > 1
            gbp = [gbp bpi(end)];
            bpd = [bpd 1];
        else
            gbp = [gbp bpi(end)-count+1];
            bpd = [bpd count];
        end
    else
        gbp = bpi;
    end
bpi = gbp;
else
    bpi = [];
    bpd = [];
end
