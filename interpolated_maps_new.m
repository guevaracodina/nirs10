function [H thz] = interpolated_maps_new(Z,W,C,Q,xCon,f1,erdf,hb,H)
stat_map = C.stat_map;
index_mask = Q.index_mask;
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
F.hb = hb;
if f1 > 0
    strf = ['_S' int2str(f1)];
    strf_fig = [' S' int2str(f1)];
else
    strf = '';
    strf_fig = '';
end
filestr = [num2str(p_value) '_' spec_hemi '_' hb strf];
filestr_fig = [num2str(p_value) ' ' spec_hemi ' ' hb strf_fig];

%CF: copy figure structure
CF.GInv = GInv;
CF.split = split;
CF.nC = nC;
if ~W.Avg
    LKC = [Q.L0 Q.L1 Q.L2];
else
    Z.output_unc = 1;
    LKC = [1 1 1];
end
%loop over contrasts
try
    nCl = 0;
    thz = cell(1,nC);
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
        
        s_map = zeros(s1, s2);
        tstr = xCon(c1).STAT;
        F.tstr = tstr;
        s_map(index_mask) = squeeze(stat_map(c1,index_mask));
        F.s_map = s_map;
        F.eidf = xCon(c1).eidf;
        %Positive responses
        %names
        F.contrast_info = [filestr '_Pos' xCon(c1).name];
        F.contrast_info_for_fig = [filestr_fig ' Pos' xCon(c1).name];
        F.contrast_info_both = [filestr xCon(c1).name]; %same for Pos and Neg, used for combined figures
        F.contrast_info_both_for_fig = [filestr_fig xCon(c1).name]; %same for Pos and Neg, used for combined figures
        CF.Z = Z; %horrible, but needed for nirs_copy_figure (near line 17)
        %uncorrected - if we generate inverted responses, or HbO or HbT
        %**************************************************************
        %Add DF to output to save the threshold value
        %Ke Peng, 2012-11-15
        %**************************************************************
        if GInv || strcmp(hb,'HbO') || strcmp(hb,'HbT')
            if Z.output_unc
                DF = nirs_draw_figure(2,F,W,Z,LKC); %DF: draw figure structure: handles to figure, axes and colorbar
                %copy figure structure
                if GFIS
                    H{nCl} = nirs_copy_figure(H{nCl},DF,CF,c1,hb,1,tstr,0,Z.write_neg_pos);
                end
            end
            if ~W.Avg
                DF = nirs_draw_figure(3,F,W,Z,LKC);
                if GFIS
                    H{nCl} = nirs_copy_figure(H{nCl},DF,CF,c1,hb,1,tstr,1,Z.write_neg_pos);
                end
                if isfield(DF, 'th_z')
                    thz{c1}.positive_thz = DF.th_z;
                end
            end
        end
        if  tstr == 'T' %for F-stat, do not invert the map - always positive
            %repeat for negative contrasts
            F.s_map = -s_map; %only used for generating Neg contrast figures
        end
        
        F.contrast_info = [filestr '_Neg' xCon(c1).name];
        F.contrast_info_for_fig = [filestr_fig ' Neg' xCon(c1).name];
        if GInv || strcmp(hb,'HbR')
            if Z.output_unc
                DF = nirs_draw_figure(2,F,W,Z,LKC);
                if GFIS, H{nCl} = nirs_copy_figure(H{nCl},DF,CF,c1,hb,0,tstr,0,Z.write_neg_pos); end
            end
            if ~W.Avg
                DF = nirs_draw_figure(3,F,W,Z,LKC);
                if GFIS, H{nCl} = nirs_copy_figure(H{nCl},DF,CF,c1,hb,0,tstr,1,Z.write_neg_pos); end
                if isfield(DF, 'th_z')
                    thz{c1}.negative_thz = DF.th_z;
                end
            end
        end
    end
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
    disp('Problem in interpolated_maps')
end