function [pE,pC] = nirs_set_priors(priorFile,Model,m)
pE = [];
pC = [];
if ~isempty(priorFile)
    dirnow = pwd;
    [pdir1 pfil1 pext1] = fileparts(priorFile);
    cd(pdir1);
    eval(['[pE,pC] = ' pfil1 '(m,3);']);
    cd(dirnow);
end
if isempty(pE) || isempty(pC)
    switch Model
        case 0 %Buxton-Friston
            %[pE,pC] = spm_hdm_priors_YO(m);%PP pour Michèle.
            [pE,pC] = spm_hdm_priors(m,3);
        case 1 %Zheng-Mayhew
            [pE,pC] = nirs_hdm_priors_ZM(m,5);   %MODIFIER
        case 2 %Huppert1
            [pE,pC] = nirs_hdm_priors_Huppert1(m,3);
        otherwise
            [pE,pC] = spm_hdm_priors(m,3);
    end
end