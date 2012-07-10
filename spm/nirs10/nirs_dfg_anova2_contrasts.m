function anova2_contrasts = nirs_dfg_anova2_contrasts(def_con)
anova2_contrasts         = cfg_entry;
anova2_contrasts.name    = 'Contrasts to include';
anova2_contrasts.tag     = 'anova2_contrasts';
anova2_contrasts.strtype = 'r';
anova2_contrasts.num     = [1 Inf];
anova2_contrasts.val     = {[1:def_con]};
anova2_contrasts.help    = {'Specify which contrasts to include in the anova'}';
