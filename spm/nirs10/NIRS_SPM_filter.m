function [NIRS d bpi bpd] = NIRS_SPM_filter(NIRS,K,d,DSon,...
                            df,bpi,bpd,markers_available)
%inputs:
%NIRS structure
%K: filter structure from NIRS_SPM
%d: data, represented as a channels X time matrix
%DSon: Boolean, if downsampling
%df: downsampling factor
%bpi: indices of bad markers, in data points
%bpd: durations of bad markers, in data points
%markers_available: Boolean
group_consecutive_markers = 0; %boolean
if markers_available
    %loop over segments and filter separately
    %group consecutive markers, as an option
    if group_consecutive_markers
        [bpi bpd] = find_marker_durations(bpi);
    end
    for i=1:length(bpi)+1
        try
            switch i
                case 1
                    %interval before first marker
                    e=d(:,1:bpi(1)-1);
                case length(bpi)+1
                    %interval after last marker
                    e=d(:,bpi(end)+bpd(end):end);
                otherwise
                    %all other intervals
                    e=d(:,bpi(i-1)+bpd(i-1):bpi(i)-1);
            end
        catch
            e = [];
        end
        
    end %end for i
    %Filter
    try
        e = spm_filter_HPF_LPF_WMDL(K,e')'; %don't forget to transpose!
    catch
        e = [];
    end
    %reconstitute d
    try
        switch i
            case 1
                %interval before first marker
                d(:,1:bpi(1)-1) = e;
            case length(bpi)+1
                %interval after last marker
                d(:,bpi(end)+bpd(end):end) = e;
            otherwise
                %all other intervals
                d(:,bpi(i-1)+bpd(i-1):bpi(i)-1) = e;
        end
    catch
        %do nothing
        disp(['Unable to filter segment ' int2str(i)]);
    end
else
    %no markers - filter over all data
    y = spm_filter_HPF_LPF_WMDL(K,d');
    d = y';
end
if DSon
    d = d(:,1:df:end);
    %how to downsample the markers?
    if markers_available
        nbpi = [];
        %loop over markers
        for i=1:length(bpi)
            %marker position after downsampling
            k = rem(bpi(i),df);
            %add to new marker list
            nbpi = [nbpi k];
            %if marker duration is more than df or even a multiple,
            %add one marker for each multiple
            for j=1:length(rem(bpd(i),df))
                k = k+1;
                nbpi = [nbpi k];
            end
        end
        %replace old markers by new, after grouping consecutive markers
        [bpi bpd] = find_marker_durations(nbpi);
    end
    %will need to change notation here - this is for
    %NIRS_SPM, which uses nirs_data for concentrations
    NIRS.nirs_data.cL.type = 'Gaussian';
    NIRS.nirs_data.K = K;
end