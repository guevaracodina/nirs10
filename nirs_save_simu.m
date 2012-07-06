function S = nirs_save_simu(M,it)
%for simulations: store results in place of subjects
if M.S.simuOn
    S.EpS = M.P;
    S.pE = M.pE;
    S.Ep = M.Ep;
    S.Cp = M.Cp;
    S.F = M.F;
    if it == 1
        S.pA = M.S.pA;
    end
end