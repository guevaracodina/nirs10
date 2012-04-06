function [TOPO H A] = call_figure_2anova(TOPO,H,Z,W,F,A,CF,v1,h1,hb,Astr,z1,LKC)
try
    %Greenhouse-Geisser and Hyunh-Feldt corrections
    if isfield(A,'Eps')
        %Recommendation: If epsilon is >0.75, use the Huynh-Feldt correction.
        %If epsilon is <0.75, or nothing is known about sphericity at all,
        %use the Greenhouse-Geisser correction
        switch z1
            %find effective position in list of epsilons
            case 1 %'intAB'
                ze = 3;
            case {2,4} %'mainA'
                ze = 1;
            case {3,5} %'mainB'
                ze = 2;
        end
        switch Z.CorrectionMethod
            case 0
                %use Greenhouse-Geisser
                epscor = A.Eps.EpsGG(ze);
            case 1
                if A.Eps.EpsGG(ze) >= 0.75
                    %use Huynh-Feldt
                    epscor = A.Eps.EpsHF(ze);
                else
                    %use Greenhouse-Geisser
                    epscor = A.Eps.EpsGG(ze);
                end
            case 2 %none
                epscor = 1;
        end
        %reduction in the number of degrees of freedom
        F.erdf = epscor*A.erdf;
        F.eidf = epscor*A.eidf;
        Astr = [Astr '_eps_' num2str(epscor,'%0.2f')];
        A.epscor = epscor;
    else
        F.erdf = A.erdf;
        F.eidf = A.eidf;
    end
    A.strA = Astr;
    F.tstr = 'F'; %tstr;
    F.hb = hb;
    filestr = [num2str(Z.p_value) '_' W.spec_hemi '_' hb '_' Astr];
    filestr_fig = [num2str(Z.p_value) ' ' W.spec_hemi ' ' hb ' ' Astr];
    info1 = filestr;
    info_for_fig1 = filestr_fig;
    F.contrast_info = info1;
    F.contrast_info_for_fig = info_for_fig1;
    F.contrast_info_both = filestr; %same for Pos and Neg, used for combined figures
    F.contrast_info_both_for_fig = filestr_fig; %same for Pos and Neg, used for combined figures
    F.s_map = A.Tmap;
    CF.contrast_info_for_fig = F.contrast_info_for_fig;
    CF.contrast_info_both_for_fig = F.contrast_info_both_for_fig;
    
    if Z.output_unc
        DF = nirs_draw_figure(7,F,W,Z,[]);
        H = nirs_copy_figure(H,DF,CF,1,hb,1,F.tstr,0,0);
    end
    if ~isempty(DF)  
        A.th_z_unc = DF.th_z;
    end
    if Z.LKC
        DF = nirs_draw_figure(5,F,W,Z,LKC);
        H = nirs_copy_figure(H,DF,CF,1,hb,1,F.tstr,1,0);
    end
    if ~isempty(DF)
        A.th_z_cor = DF.th_z;
    end
    if z1 <= 3
        TOPO.v{v1}.group.hb{h1}.c{2*z1-1} = A;
    else
        if z1 == 4
            TOPO.v{v1}.group.hb{h1}.effAonB{2*A.y1-1} = A;
        else
            if z1 == 5
                TOPO.v{v1}.group.hb{h1}.effBonA{2*A.y1-1} = A;
            end
        end
    end
catch exception2
    disp(exception2.identifier);
    disp(exception2.stack(1));
end