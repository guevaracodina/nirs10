function outputdatafigures = nirs_dfg_outputdatafigures
outputdatafigures      = cfg_menu;
outputdatafigures.tag  = 'outputdatafigures';
outputdatafigures.name = 'Output data figures';
outputdatafigures.labels = {'Yes', 'No'};
outputdatafigures.values = {1,0};
outputdatafigures.val{1}  = 1;
outputdatafigures.help = {'Output data figures.' }';
