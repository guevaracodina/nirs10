function [xCon nC Z] = nirs_fix_for_AvgInterpBetaMode_2anova(nC,xCon,big_TOPO,v1,Z,NIRS)
if nC == 0 || isempty(xCon{1})
    clear xCon
    Z.AvgInterpBetaMode = 1;
    %quick fix for case Z.AvgInterpBetaMode (akin to variable used in liom_group)
    nC = 1;
    Ns = length(big_TOPO{1}.v{v1}.s); %number of sessions
    %need to factorize sessions into Z.anova2_sessions x Z.anova2_contrasts
    effNs = length(Z.anova2_sessions);
    effNc = length(Z.anova2_contrasts);
    s1ctr = 0;
    if Ns == effNs * effNc
        for s1 = 1:effNs
            for c1 = 1:effNc;
                s1ctr = s1ctr + 1;
                try
                    load(NIRS.SPM{1});
                    xCon{s1}(c1).name = SPM.xXn{s1ctr}.Sname;
                catch
                    xCon{s1}(c1).name = ['A' gen_num_str(s1ctr,2)]; %could put the user-specified names
                end
                xCon{s1}(c1).STAT = 'T';
                xCon{s1}(c1).SessionNumber = s1; %check
                xCon{s1}(c1).eidf = 1;
            end
        end
    end
end