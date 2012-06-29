function channel_choice = nirs_dfg_channel_choice
all_channels         = cfg_branch;
all_channels.tag     = 'all_channels';
all_channels.name    = 'All channels';
all_channels.val     = {};
all_channels.help    = {'All channels will be processed'};

selected_channels      = cfg_entry;
selected_channels.tag  = 'selected_channels';
selected_channels.name = 'Enter list of channels';
selected_channels.strtype  = 'r';
selected_channels.num = [1 Inf];
selected_channels.val{1} = 1;
selected_channels.help = {'Enter list of channels to process.'};

select_channels         = cfg_branch;
select_channels.tag     = 'select_channels';
select_channels.name    = 'Select channels';
select_channels.val     = {selected_channels};
select_channels.help    = {'Choose some channels to be processed'};

channel_choice        = cfg_choice;
channel_choice.name   = 'Choose channel selection method';
channel_choice.tag    = 'channel_choice';
channel_choice.values = {all_channels,select_channels};
channel_choice.val    = {all_channels};
channel_choice.help   = {'Choose channel selection method'}';