function [useSeg XYZ_mm_f] = nirs_cleanup_c1(fSeg,skinOn)
if spm_existfile(fSeg)
    try
        V = spm_vol(fSeg);
        [Y0 XYZ] = spm_read_vols(V);
        %clean up the segmented file
        if skinOn
            layer = 5;
        else
            %cortex
            layer = 1;
        end
        XYZ_mm_f = XYZ(:,Y0==layer); %select c1 layer or skin layer
        Nv = size(XYZ_mm_f,2);
        vxXYZ_mm_f = V.mat\[XYZ_mm_f; ones(1,Nv)];
        A = nirs_clusters(round(vxXYZ_mm_f(1:3,:)),6);
        nC = max(A);
        cs = zeros(1,nC);
        for c1=1:nC
            cs(c1) = sum(A == c1);
        end
        [csort cidx] = sort(cs,'descend');
        min_size = 1e4;
        idx_keep = cidx(csort>min_size);
        %remove those voxels from the c1 layer
        Ac = zeros(1,Nv);
        for c1=1:length(idx_keep)
            Ac = Ac + (A == idx_keep(c1));
        end
        XYZ_mm_f = XYZ_mm_f(:,logical(Ac));
        useSeg = 1;
    catch 
        useSeg = 0;
        XYZ_mm_f = [];
    end
end