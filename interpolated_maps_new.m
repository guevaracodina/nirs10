function H = interpolated_maps_new(Z,W,C,Q,xCon,f1,erdf,hb,H)
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
LKC = [Q.L0 Q.L1 Q.L2];
%loop over contrasts
try
    for c1=1:nC
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
        
        %uncorrected - if we generate inverted responses, or HbO or HbT
        if GInv || strcmp(hb,'HbO') || strcmp(hb,'HbT')
            if Z.output_unc
                DF = nirs_draw_figure(2,F,W,Z,LKC); %DF: draw figure structure: handles to figure, axes and colorbar
                %copy figure structure
                if GFIS, H = nirs_copy_figure(H,DF,CF,c1,hb,1,tstr,0,Z.write_neg_pos); end
            end
            DF = nirs_draw_figure(3,F,W,Z,LKC);
            if GFIS, H = nirs_copy_figure(H,DF,CF,c1,hb,1,tstr,1,Z.write_neg_pos); end
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
                if GFIS, H = nirs_copy_figure(H,DF,CF,c1,hb,0,tstr,0,Z.write_neg_pos); end
            end
            DF = nirs_draw_figure(3,F,W,Z,LKC);
            if GFIS, H = nirs_copy_figure(H,DF,CF,c1,hb,0,tstr,1,Z.write_neg_pos); end
        end
    end
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
    disp('Problem in interpolated_maps')
end