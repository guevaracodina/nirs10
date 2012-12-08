function Pp_c1_rmm = nirs_coreg_optodes_unnormalized(Pp_rmm,Pvoid,...
    p_cutoff,fc1,fSeg,coreg_ch,dir_coreg)

% This function is developped specially for coregistering optodes positions
% from scalp to c1 grey matter layer. The coordinates of optodes should
%==========================================================================
% Inputs:
% Pp_mm:  optodes coordinates in mm;
% Pvoid: matrix indicating functional optodes and void optodes;
% Outputs:
% Pp_c1_mm: unnormalised optodes coordinates in mm (matched to a point on c1)
%==========================================================================
% Molecular and Optical Imaging Laboratory
% Ecole Polytechnique de Montreal
% Ke Peng
% 2012-09-21
useSeg = 0;
if ~isempty(fSeg)
    [useSeg XYZ_mm_f] = nirs_cleanup_c1(fSeg);
end
if ~useSeg
    %Extract c1 image
    if spm_existfile(fc1)
        V = spm_vol(fc1);
        %[dir1 fil1] = fileparts(fc1);
    else
        disp(['Cannot find c1 image for optode re-coregistration at location: ' fc1]);
        Pp_c1_rmm = [];
        return;
    end
    
    Pp_c1_rmm = zeros(size(Pp_rmm));
    [Y0 XYZ] = spm_read_vols(V);
    
    XYZ_mm_f = XYZ;
    XYZ_mm_f(:,Y0<=p_cutoff) = [];
end
dis = [];
d = pdist2(XYZ_mm_f',Pp_rmm');

for i = 1 : size(Pp_rmm,2)
    if coreg_ch || ~Pvoid(1,i)
        use_sort = 0;
        if use_sort
            [sd Idx] = sort(d(:,i));
            xsmallest = 1;
            a = sd(xsmallest); %xsmallest, to avoid "rebel" voxels
            b = Idx(xsmallest);
        else
            [a b] = min(d(:,i));
        end
        dis = [dis d(b,i)];
        Pp_c1_rmm(1:3,i) = XYZ_mm_f(:,b);
        Pp_c1_rmm(4,i) = 1;
    else
        Pp_c1_rmm(1:3,i) = 0;
        dis = [dis -1];
    end
end
if coreg_ch
    str = '_channels';
else
    str = '_optodes';
end
save(fullfile(dir_coreg,['coreg_c1_distance' str '.mat']),'dis');



