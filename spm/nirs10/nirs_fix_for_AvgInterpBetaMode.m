function [xCon nC Z] = nirs_fix_for_AvgInterpBetaMode(nC,xCon,big_TOPO,v1,Z,NIRS)
if nC == 0 || isempty(xCon{1})
    clear xCon
    Z.AvgInterpBetaMode = 1;
    %quick fix for case Z.AvgInterpBetaMode (akin to variable used in liom_group)
    nC = 1;
    Ns = length(big_TOPO{1}.v{v1}.s); %number of sessions
    for s1 = 1:Ns
        try
            load(NIRS.SPM{1});
            xCon{s1}.name = SPM.xXn{s1}.Sname;
        catch
            xCon{s1}.name = ['A' gen_num_str(s1,2)]; %could put the user-specified names
        end
        xCon{s1}.STAT = 'T';
        xCon{s1}.SessionNumber = s1;
        xCon{s1}.eidf = 1;
    end
end