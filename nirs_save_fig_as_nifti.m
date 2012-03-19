function nirs_save_fig_as_nifti(Z,W,pathn,combinedfig,tstr,str_cor,contrast_info,F)
%Save as nifti
if Z.save_nifti_contrasts && ( Z.write_neg_pos || combinedfig || tstr == 'F' )
    pathnii = fullfile(pathn,'nii');
    if ~exist(pathnii,'dir'),mkdir(pathnii); end
    %NP = not permuted
    filen5 = fullfile(pathnii,[tstr '_' str_cor '_' contrast_info 'NP.nii']);
    filen3 = fullfile(pathnii,[tstr '_' str_cor '_' contrast_info '.nii']);
    %note it is the contrast that should be written, not the T or F-stat maps
    M = [[0 1;-1 0] zeros(2); zeros(2) eye(2)];
    if strcmp(tstr,'T')
        V = nirs_create_vol(filen5,...
            [W.s1 W.s2 1], [16,0], [1;0;352],M, F.con);
    else
        V = nirs_create_vol(filen5,...
            [W.s1 W.s2 1], [16,0], [1;0;352],M, F.ess);
    end
    %             end
    % test : on se trompe pas pour le milieu
    % 1:ventral, 2:dorsal, 3:right, 4:left, 5:frontal, 6:occipital
    switch W.side_hemi
        case 2% s2 x vx and s1 z vx
            dim = [W.s2 1 W.s1];
            M = [[0 0 1;1 0 0;0 1 0;0 0 0] zeros(4,1)];
            vecta = -M(1:3,1:3)*[round(dim(1)/2);1;round(dim(3)/2)];
            Fcon = permute(F.con,[2,3,1]);
            Fess = permute(F.ess,[2,3,1]);
            
        case 3% s2 -x vx and s1 y vx
            dim = [W.s2 W.s1 1];
            M = [[0 0 1;1 0 0;0 1 0;0 0 0] zeros(4,1)];
            vecta = -M(1:3,1:3)*[round(dim(1)/2);round(dim(2)/2);1];
            Fcon = permute(F.con,[2,1,3]);
            Fess = permute(F.ess,[2,1,3]);
        case 4% s2 x vx and s1 y vx
            dim = [W.s2 W.s1 1];
            M = [[0 0 1;-1 0 0;0 1 0;0 0 0] zeros(4,1)];%%% matrice de base
            vecta = -M(1:3,1:3)*[round(dim(1)/2);round(dim(2)/2);1];
            Fcon = permute(F.con,[2,1,3]);
            Fess = permute(F.ess,[2,1,3]);
            
        case 5% s2 -z vx and s1 y vx
            dim = [1 W.s1 W.s2];
            M = [[0 0 1;-1 0 0;0 -1 0;0 0 0] zeros(4,1)];
            vecta = -M(1:3,1:3)*[1;round(dim(1)/2);round(dim(2)/2)];
            Fcon = permute(F.con,[3,1,2]);
            Fess = permute(F.ess,[3,1,2]);
            %case 6% s2 z vx and s1 y vx
        otherwise % cas sans interet
    end
    M(:,4) = [vecta;1];
    % test : use the V.mat of the T1
    %         path_T1 = [fileparts(fileparts(pathn)) '\T1'];
    %         V = spm_vol([path_T1 '\T1.nii']);
    %         M(1:3,1:3) = V.mat(1:3,1:3)./norm(V.mat(1:3,1:3));
    if strcmp(tstr,'T')
        V = nirs_create_vol(filen3,...
            dim, [16,0], [1;0;352],M, Fcon);
        %                 [W.s1 W.s2 1], [16,0], [1;0;352],M, F.con);
        
    else
        V = nirs_create_vol(filen3,...
            dim, [16,0], [1;0;352],M, Fess);
        %                 [W.s1 W.s2 1], [16,0], [1;0;352],M, F.ess);
        
    end
end