function Pp_c1_wmm = nirs_coreg_optodes(Pp_wmm,Pvoid,coreg_optodes,T1_loc)

% This function is developped specially for coregistering optodes positions
% from scalp to c1 grey matter layer. The coordinates of optodes should be
% all normalised.
%==========================================================================
% Inputs:
    % Pp_wmm: Normalised optodes coordinates in mm;
    % Pvoid: matrix indicating functional optodes and void optodes;
    % coreg_optodes: configurations from previous batch;
    % image_loc: string indicating the location of c1 image;
% Outputs:
    % Pp_c1_wmm: Normalised optodes coordinates in mm (matched to a point on c1)
%==========================================================================
% Molecular and Optical Imaging Laboratory
% Ecole Polytechnique de Montreal
% Ke Peng
% 2012-09-21

if ~isstruct(coreg_optodes)
    disp('ERROR! Config for coregistering optodes onto c1 is not a struct!');
    Pp_c1_wmm = [];
    return;
end

[dir0 fil0] = fileparts(T1_loc);
[files,dirs] = spm_select('FPList',dir0,'wc1*');%First try to find wc1 file
found = [];
for f0=1:size(files,1)
    t1 = files(f0,:);
    [dir1 fil1 ext1] = fileparts(t1);
    if strcmp(deblank(ext1),'.nii')
        %should be a good c1 file
        found = t1;
        break;
    end
end

if ~isempty(found)
    image_loc = deblank(found);
else
    disp(['Cannot find wc1 image at location' dir1]);
    Pp_c1_wmm = [];
    return;
end

%Extract c1 image
if spm_existfile(image_loc)
    V = spm_vol(image_loc);
else
    disp(['Cannot find c1 image for optode re-coregistration at location: ' image_loc]);
    Pp_c1_wmm = [];
    return;
end


if isfield(coreg_optodes.coreg_c1_method, 'coreg_c1_ball')

    r_op = coreg_optodes.coreg_c1_method.coreg_c1_ball.coreg_radius;
    r_op_retry = coreg_optodes.coreg_c1_method.coreg_c1_ball.coreg_radius_retry;
    p_cutoff = coreg_optodes.coreg_c1_method.coreg_c1_ball.coreg_ball_cutoff;



    %Do matching (optodes position on scalp => nearest point on c1)
    Pp_c1_wmm = zeros(size(Pp_wmm));

    for i = 1 : size(Pp_wmm,2)
        if ~Pvoid(1,i)                     %Do search for all the optodes
        %including the void ones
            last_radius = r_op;
            found_op = search_coreg_point_c1(Pp_wmm(:,i),last_radius,0,V);

            while isempty(found_op)
                last_radius = last_radius + r_op_retry;
                disp(['Last iteration did not find any corresponding point for optode No.' int2str(i)]);
                disp('Start another iteration...')
    %             method1 = 0;
    %             if method1
                found_op = search_coreg_point_c1(Pp_wmm(:,i),last_radius,(last_radius-r_op_retry),V,p_cutoff);
    %             else
    %                 [Y,XYZ] = spm_read_vols(V);
    %                 tmpVo = V.mat\[Pp_wmm(:,i);1];
    %                 d = (XYZ-tmpVo).^2
    %                 [a b] = min(d);
    %                 nearestPoint = XYZ(b);
    %                 
    %             end
    %             
            end
            Pp_c1_wmm(1:3,i) = found_op;
            Pp_c1_wmm(4,i) = 1;
            disp(['Coregistration for optode No.' int2str(i) ' finished.']);
        else
            Pp_c1_wmm(:,i) = Pp_wmm(:,i);
        end
    end
else
    p_cutoff = coreg_optodes.coreg_c1_method.coreg_c1_cutoff.coreg_global_cutoff;
    
    Pp_c1_wmm = zeros(size(Pp_wmm));
    [Y0 XYZ] = spm_read_vols(V);
    
    XYZ_mm_filtered = [];
    xt = size(Y0,1);
    yt = size(Y0,2);
    zt = size(Y0,3);
    for i0 = 1:xt
        for y0 = 1:yt
            for z0 = 1:zt
                if Y0(i0,y0,z0) > p_cutoff
                    XYZ_mm_filtered = [XYZ_mm_filtered XYZ(:,i0+(y0-1)*xt+(z0-1)*xt*yt)];
                end
            end
        end
    end
    
    %XYZ_mm_filtered = V.mat*[XYZ_filtered;ones(1,size(XYZ_filtered,2))];
    dis = [];
    for i = 1 : size(Pp_wmm,2)
        
        if ~Pvoid(1,i)
     
            %tmpVo = V.mat\[Pp_wmm(:,i);1];
            d0(1,:) = XYZ_mm_filtered(1,:) - Pp_wmm(1,i);
            d0(2,:) = XYZ_mm_filtered(2,:) - Pp_wmm(2,i);
            d0(3,:) = XYZ_mm_filtered(3,:) - Pp_wmm(3,i);
            d0 = d0.^2;
            d = d0(1,:)+d0(2,:)+d0(3,:);
            [a b] = min(d);
            dis = [dis sqrt(d(b))];
            Pp_c1_wmm(1:3,i) = XYZ_mm_filtered(1:3,b);
            Pp_c1_wmm(4,i) = 1;
            %disp(['Coregistration for optode No.' int2str(i) ' finished.']);
        else
            Pp_c1_wmm(:,i) = Pp_wmm(:,i);
            dis = [dis -1];
        end
    end
    
    save([dir1 '\coreg_c1_distance.mat'],'-mat','dis');
  
end


function P_c1_wmm = search_coreg_point_c1(P_wmm,l_r,p_r,V,p_cutoff)

P_tmp_wmm = [];
P_tmp_wvx = [];

P_all_contact = [];

for x_loop = -l_r : l_r
    
    for y_loop = -l_r:l_r
        
        for z_loop = -l_r:l_r
            
            if x_loop^2 + x_loop^2 + x_loop^2 <= l_r*l_r %Ball assumption
                
                if x_loop^2 + x_loop^2 + x_loop^2 <= p_r*p_r %If this coordinate has been checked in previous loops, then jump to next coordinates
                    continue;
                else
                    x0 = P_wmm(1,1) + x_loop;
                    y0 = P_wmm(2,1) + y_loop;
                    z0 = P_wmm(3,1) + z_loop;
                    
                    P_tmp_wmm = [x0;y0;z0;1];

                    %Change optode coordinates mm => voxels
                    P_tmp_wvx = V.mat\P_tmp_wmm;

                    inten_P = spm_sample_vol(V,P_tmp_wvx(1,1),P_tmp_wvx(2,1),P_tmp_wvx(3,1),1);

                    if abs(inten_P) > p_cutoff
                        d0 = sqrt((x_loop)^2 + (y_loop)^2 + (z_loop)^2);%LS method
                        P_all_contact = [P_all_contact [x0;y0;z0;d0]];               
                    end
                end
            end
        end
    end
end

if isempty(P_all_contact)
    P_c1_wmm = [];
else
    [i0 j0] = find(P_all_contact(4,:) == min(P_all_contact(4,:)));
    P_c1_wmm = [P_all_contact(1,j0(1));P_all_contact(2,j0(1));P_all_contact(3,j0(1))]; %If more than one corresponding point, simly take the first one.
end

