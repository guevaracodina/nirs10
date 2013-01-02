function cortex_projection_method = nirs_dfg_cortex_projection_method

project_liom         = cfg_branch;
project_liom.tag     = 'project_liom';
project_liom.name    = 'Project onto cortex using liom1 method (PP)';
project_liom.val     = {};
project_liom.help    = {'Project onto cortex using liom1 method.'};

project_liom_Ke         = cfg_branch;
project_liom_Ke.tag     = 'project_liom_Ke';
project_liom_Ke.name    = 'Project onto cortex using liom2 method (Ke)';
project_liom_Ke.val     = {};
project_liom_Ke.help    = {'Project onto cortex using liom2 method.'};

project_Korean         = cfg_branch;
project_Korean.tag     = 'project_Korean';
project_Korean.name    = 'Project onto cortex using Korean method';
project_Korean.val     = {};
project_Korean.help    = {'roject onto cortex using Korean method.'
    'This calculation is done using a specific Korean template.'}';

cortex_projection_method        = cfg_choice;
cortex_projection_method.name   = 'Cortex projection choice';
cortex_projection_method.tag    = 'cortex_projection_method';
cortex_projection_method.values = {project_liom,project_liom_Ke,project_Korean};
cortex_projection_method.val    = {project_Korean};
cortex_projection_method.help   = {'Projection onto cortex.'}';
