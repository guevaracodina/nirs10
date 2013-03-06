function TOPO = fill_TOPO(TOPO,C,side_hemi,f1,hb,thz)
%Fills TOPO structure with interpolated beta, cov_beta, and F-contrast
hbi = inv_get_chromophore(hb);
%CONVENTION:
%'g': group of sessions
%'group': group analysis over subjects (see liom_group)
%'s': individual sessions
if f1 == 0
    fg = 'g';
    f1 = 1;
else
    fg = 's';
end
if isfield(C,'sum_kappa')
    TOPO.v{side_hemi}.(fg){f1}.hb{hbi}.sum_kappa = C.sum_kappa;
    TOPO.v{side_hemi}.(fg){f1}.hb{hbi}.c_interp_beta = C.c_interp_beta;
    TOPO.v{side_hemi}.(fg){f1}.hb{hbi}.c_cov_interp_beta = C.c_cov_interp_beta;
    TOPO.v{side_hemi}.(fg){f1}.hb{hbi}.c_interp_F = C.c_interp_F;
    TOPO.v{side_hemi}.(fg){f1}.hb{hbi}.c_interp_T = C.c_interp_T;
    TOPO.v{side_hemi}.(fg){f1}.hb{hbi}.c_interp_ess = C.c_interp_ess;
    TOPO.v{side_hemi}.(fg){f1}.hb{hbi}.c_interp_ess0 = C.c_interp_ess0;
else
    try %skips the case Z.DoStats == 0
    TOPO.v{side_hemi}.(fg){f1}.hb{hbi}.stat_map = C.stat_map;
    end
    TOPO.v{side_hemi}.(fg){f1}.hb{hbi}.beta_map = C.beta_map;
    %TOPO.v{side_hemi}.(fg){f1}.hb{hbi}.cov_map = C.cov_map;
end
if ~isempty(thz)
    for h0 = 1:length(thz)
        if isfield(thz{h0},'positive_thz') || isfield(thz{h0},'negative_thz')
            if isfield(TOPO.v{side_hemi}.(fg){f1}.hb{hbi},'th_z')
                h1 = length(TOPO.v{side_hemi}.(fg){f1}.hb{hbi}.th_z);
            else
                h1 = 0;
            end
            TOPO.v{side_hemi}.(fg){f1}.hb{hbi}.th_z{h1+1} = thz{h0};
            %Record all the H that contain a "DF.thz" field
        end
    end
end
TOPO.v{side_hemi}.(fg){f1}.hb{hbi}.hb = hb;