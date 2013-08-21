function [p_value t_value flag first_k_value stop_k_value] = nirs_2D_Bonferroni_FDR_threshold(CType,tmap,rmap,ch,q_value,erdf,OP)
% This function is used to calculate the p-threshold and the t-threshold
% for a 2D projected contrast map based on voxel-wise methods.
% Ke Peng, 2013-08-11

tmap = tmap .* rmap; %applying mask

switch CType
    case 1 % Pixel-wise Bonferroni correction
        N = length(find(rmap));
        p_value = q_value/N;
        t_value = spm_invTcdf(p_value,erdf);
        flag = 1;
        stop_k_value = -1; %Not used
        first_k_value = N;
        
    case 2 % Pixel-wise FDR correction
        FDR_type = OP.fdrtype;
        t_type = OP.t_type;
        
        z0 = tmap(find(tmap));
        Z = z0';
        switch t_type
            case 1 %One tailed
                p_pixel = spm_Tcdf(Z,erdf);
            case 2 %Two tailed
                p_pixel = spm_Tcdf(Z,erdf)*2;
            otherwise
        end
        for zi = 1 : length(Z)
            if Z(zi) > 0
                p_pixel(zi) = 1 - p_pixel(zi); %for positive t-value, use 1-p_value
            end
        end
        [psort pidx] = sort(p_pixel);
        N = length(psort); %Number of pixels inside the boundray
        %Perform FDR
        [p_value flag first_k_value stop_k_value] = nirs_get_FDR_threshold(FDR_type,q_value,N,psort);
        z_idx = pidx(stop_k_value);
        t_value = Z(z_idx);                
%                 t_value = spm_invTcdf(p_value/2,erdf);

    case 3 % Peak-FDR correction
        FDR_type = OP.fdrtype;
        u_thz = OP.u_thz;
        t_type = OP.ttype;
        
        switch t_type
            case 1  % One-sided t-test
                if ch == 2 % HbR
                    l0 = find(tmap > -u_thz);
                    tmap2 = tmap;
                    tmap2(l0) = 0; %Apply the first threshold
                    if max(max(abs(tmap2))) == 0
                        disp('2D Peak-FDR: No peak passes the initial threshold.');
                        flag = 0;
                        p_value = 0;t_value = 0;
                        first_k_value = -1;stop_k_value = -1;
                        return
                    end                        
                    [Zt Mt] = find_peak_2D(abs(tmap2)); %All t-values to be positive (for peak localization)
                    Zt = -Zt; % change back to negative
                    p_peak = spm_Tcdf(Zt,erdf);
                else % HbO or HbT
                    l0 = find(tmap < u_thz);
                    tmap2 = tmap;
                    tmap2(l0) = 0; %Apply the first threshold
                    if max(max(abs(tmap2))) == 0
                        disp('2D Peak-FDR: No peak passes the initial threshold.');
                        flag = 0;
                        p_value = 0;t_value = 0;
                        first_k_value = -1;stop_k_value = -1;
                        return
                    end  
                    [Zt Mt] = find_peak_2D(tmap2);
                    p_peak = 1 - spm_Tcdf(Zt,erdf); % for positive t-value, use 1-p_value
                end
            case 2  % Two-sided t-test
                l0 = find(abs(tmap) < u_thz);
                tmap2 = tmap;
                tmap2(l0) = 0; %Apply the first threshold
                if max(max(abs(tmap2))) == 0
                    disp('2D Peak-FDR: No peak passes the initial threshold.');
                    flag = 0;
                    p_value = 0;t_value = 0;
                    first_k_value = -1;stop_k_value = -1;
                    return
                end                  
                [Zt Mt] = find_peak_2D(abs(tmap2)); %All t-values to be positive (for peak localization) 
                p_peak = (1-spm_Tcdf(Zt,erdf))*2; % Two-tailed p-value
%                 [psort pidx] = sort(p_peak);  %Sort in ascending order
%                 [p_value flag first_k_value stop_k_value] = nirs_get_FDR_threshold(FDR_type,q_value,length(psort),psort);
%                 t_value = spm_invTcdf(p_value/2,erdf); % Two-tailed t-value
            otherwise
        end
        [psort pidx] = sort(p_peak);  %Sort in ascending order
        %Perform FDR
        [p_value flag first_k_value stop_k_value] = nirs_get_FDR_threshold(FDR_type,q_value,length(psort),psort);
        z_idx = pidx(stop_k_value);
        t_value = Zt(z_idx);
        if flag == 1 && stop_k_value == 1
            %t_value = floor(abs(t_value)*10)/10 * (abs(t_value))/t_value; % A temporaory solution for single point problem
            t_value = t_value - 0.01*sign(t_value);
            disp('Single point passing pFDR. Imposing a less strict threshold (-0.01)');
        end
end

function [Zt Mt] = find_peak_2D(tmap)
% find peaks for a 2D projected map
[x0 y0] = find(tmap);
z0 = tmap(find(tmap));
L = [x0';y0';ones(1,length(x0'))];
Z = z0';
[Nt,Zt,Mt,At,XYZt] = spm_max(Z,L); %find peaks
