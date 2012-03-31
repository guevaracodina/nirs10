function H = interpolated_maps(Z,W,C,Q,xCon,f1,erdf,hb,H)
%This still requires clean-up; Relabel T_map for stat_map; get rid of
%c_interp_T, F, ess, ess0
sum_kappa = C.sum_kappa;
c_interp_beta = C.c_interp_beta;
c_cov_interp_beta = C.c_cov_interp_beta;
c_interp_T = C.c_interp_T;
c_interp_F = C.c_interp_F;
c_interp_ess = C.c_interp_ess;
c_interp_ess0 = C.c_interp_ess0;
s1 = W.s1;
s2 = W.s2;
spec_hemi = W.spec_hemi;
p_value = Z.p_value;
GInv = Z.GInv;
GFIS = Z.GFIS;

nC = size(xCon,2);
%Structure F
load Split
F.split = split;
F.pathn = Z.dir1;
F.erdf = erdf;
if f1 > 0
    strf = ['_S' int2str(f1)];
    strf_fig = [' S' int2str(f1)];
else
    strf = '';
    strf_fig = '';
end
filestr = [num2str(p_value) '_' spec_hemi '_' hb strf '_'];
filestr_fig = [num2str(p_value) ' ' spec_hemi ' ' hb strf_fig ' '];

%CF: copy figure structure
CF.GInv = GInv;
CF.split = split;
CF.nC = nC;

%loop over contrasts
try
    nCl = 0;
    for c1=1:nC
        if isfield(Z,'use_nCloop')
            if Z.use_nCloop               
                Z.c1eff = mod(c1,Z.nCloop);
                if Z.c1eff == 0
                    Z.c1eff = Z.nCloop;
                end
                if Z.c1eff == 1
                    nCl = nCl + 1;
                end
            end
        else
            nCl = 1;
        end
        if xCon(c1).STAT == 'T'
            impose_T_bound = 1;
            if ~isempty(c_cov_interp_beta)
                index_mask = find(squeeze(c_cov_interp_beta(c1,:,:)) ~= 0);
            else
                index_mask = Q.index_mask;
                impose_T_bound = 0;
            end
            
            s_map = zeros(s1, s2);
            if impose_T_bound
                %put a bound on low variance -- as in SPM (see ResMS bound)
                tmp = squeeze(c_cov_interp_beta(c1,index_mask));
                bound = (1e-3)^0.5*max(tmp(isfinite(tmp)));
                tmp2 = bound*ones(size(tmp));
                tmp2(tmp > bound) = tmp(tmp > bound);
                s_map(index_mask) = squeeze(c_interp_beta(c1,index_mask))./ ...
                    sqrt(tmp2); %c_interp_T(c1,index_mask);
            else
                if Z.LKC || Z.UseCorrelRes
                    s_map(index_mask) = squeeze(c_interp_T(c1,index_mask));
                else
                    s_map(index_mask) = squeeze(c_interp_beta(c1,index_mask))./ ...
                        sqrt(squeeze(c_cov_interp_beta(c1,index_mask)));
                end
            end
            tstr = 'T';
        else
            %F-stats
            if ~isempty(c_interp_F)
                index_mask = find(squeeze(c_interp_F(c1,:,:)) ~= 0);
                impose_F_bound = 1;
            else
                index_mask = Q.index_mask;
                impose_F_bound = 0;
            end
            %still use variable T_map, though it is now an F_map
            s_map = zeros(s1, s2);
            
            if impose_F_bound
                tmp = squeeze(c_interp_ess0(index_mask));%should calculate just once
                bound = 1e-3*max(tmp(isfinite(tmp)));
                tmp2 = bound*ones(size(tmp));
                tmp2(tmp > bound) = tmp(tmp > bound);
                s_map(index_mask) = squeeze(c_interp_ess(c1,index_mask))./tmp2';
            else
                s_map(index_mask) = squeeze(c_interp_F(c1,index_mask));
            end
            tstr = 'F';
        end
        %Positive responses
        %names
        F.contrast_info = [filestr '_Pos' xCon(c1).name];
        F.contrast_info_for_fig = [filestr_fig ' Pos' xCon(c1).name];
        F.contrast_info_both = [filestr xCon(c1).name]; %same for Pos and Neg, used for combined figures
        F.contrast_info_both_for_fig = [filestr_fig xCon(c1).name]; %same for Pos and Neg, used for combined figures
        
        F.s_map = s_map;
        F.eidf = xCon(c1).eidf;
        if ~isempty(sum_kappa)
            F.sum_kappa = sum_kappa(c1);
        end
        F.tstr = tstr;
        F.hb = hb;
        %contrast - what will be saved as nifti if requested
        if strcmp(tstr,'T')
            %F.con = zeros(s1,s2);
            %F.con = squeeze(c_interp_beta(c1,index_mask));
            F.con = squeeze(c_interp_beta(c1,:,:));
            F.ess = [];
        else
            F.con = [];
            %F.ess = zeros(s1,s2);
            %F.ess = squeeze(c_interp_F(c1,index_mask));
            F.ess = squeeze(c_interp_ess(c1,:,:));
        end
        if Z.LKC || Z.UseCorrelRes
            LKC = [Q.L0 Q.L1 Q.L2];
        else 
            LKC = [];
        end
        CF.Z = Z; %horrible, but needed for nirs_copy_figure (near line 17)
        %uncorrected - if we generate inverted responses, or HbO or HbT
        if GInv || strcmp(hb,'HbO') || strcmp(hb,'HbT')
            if Z.output_unc
                DF = nirs_draw_figure(2,F,W,Z,LKC); %DF: draw figure structure: handles to figure, axes and colorbar
                %copy figure structure
                if GFIS, H{nCl} = nirs_copy_figure(H{nCl},DF,CF,c1,hb,1,tstr,0,Z.write_neg_pos); end
            end
            %tube formula corrected for Tstats - or Bonferroni for Fstats
            %if tstr == 'T' %only for Tstats for now
            DF = nirs_draw_figure(3,F,W,Z,LKC);
            if GFIS, H{nCl} = nirs_copy_figure(H{nCl},DF,CF,c1,hb,1,tstr,1,Z.write_neg_pos); end
            %end
        end
        if  tstr == 'T' %for F-stat, do not invert the map - always positive
            %repeat for negative contrasts
            F.s_map = -s_map; %only used for generating Neg contrast figures
        end
        
        F.contrast_info = [filestr '_Neg' xCon(c1).name];
        F.contrast_info_for_fig = [filestr_fig ' Neg' xCon(c1).name];
        %contrast - what will be saved as nifti if requested
        if strcmp(tstr,'T')
            F.con = squeeze(c_interp_beta(c1,:,:));
            F.ess = [];
        else
            F.con = [];
            F.ess = squeeze(c_interp_ess(c1,:,:));
        end
        if GInv || strcmp(hb,'HbR')
            if Z.output_unc
                DF = nirs_draw_figure(2,F,W,Z,LKC);
                if GFIS, H{nCl} = nirs_copy_figure(H{nCl},DF,CF,c1,hb,0,tstr,0,Z.write_neg_pos); end
            end
            DF = nirs_draw_figure(3,F,W,Z,LKC);
            if GFIS, H{nCl} = nirs_copy_figure(H{nCl},DF,CF,c1,hb,0,tstr,1,Z.write_neg_pos); end
        end
    end
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
    disp('Problem in interpolated_maps')
end