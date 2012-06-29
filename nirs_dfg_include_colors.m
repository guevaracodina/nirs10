function IC = nirs_dfg_include_colors(HbO,HbR,HbT)
include_HbT      = cfg_menu;
include_HbT.tag  = 'include_HbT';
include_HbT.name = 'Include HbT';
include_HbT.labels = {'Yes','No'};
include_HbT.values = {1,0};
include_HbT.val  = {HbT};
include_HbT.help = {'Include HbT.'}';

include_HbO      = cfg_menu;
include_HbO.tag  = 'include_HbO';
include_HbO.name = 'Include HbO';
include_HbO.labels = {'Yes','No'};
include_HbO.values = {1,0};
include_HbO.val  = {HbO};
include_HbO.help = {'Include HbO.'}';

include_HbR      = cfg_menu;
include_HbR.tag  = 'include_HbR';
include_HbR.name = 'Include HbR';
include_HbR.labels = {'Yes','No'};
include_HbR.values = {1,0};
include_HbR.val  = {HbR};
include_HbR.help = {'Include HbR.'}';

IC         = cfg_branch;
IC.name     = 'Include colors';
IC.tag    = 'IC';
IC.val     = {include_HbO include_HbR include_HbT}; 
IC.help    = {'Choose colors to include.'};