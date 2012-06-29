function Model_Choice = nirs_dfg_Model_Choice
Model_Choice           = cfg_menu;
Model_Choice.name      = 'Model choice';
Model_Choice.tag       = 'Model_Choice';
Model_Choice.labels    = {'Buxton-Friston' 'Zheng-Mayhew' 'Boas-Huppert' 'BF eq. 1-2 (flow)' 'BF eq. 3-4 (bold)'};
Model_Choice.values    = {0,1,2,3,4};
Model_Choice.val       = {0};
Model_Choice.help      = {'Choose hemodynamic model: Buxton-Friston, '
    'Zheng-Mayhew, or 1-Compartment Boas-Huppert Model. '
    'The 2 last options make use of the Buxton-Friston model, but separate '
    'the linear (neural->flow, first 2 state equations) from the non-linear '
    '(flow->bold, last 2 state equations) parts.' }';
