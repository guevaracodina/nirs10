function which_session = nirs_dfg_which_session
which_session     = cfg_entry;
which_session.name    = 'Which session(s)?';
which_session.tag     = 'which_session';
which_session.strtype = 'r';
which_session.val     = {1};
which_session.num     = [0 Inf];
which_session.help    = {'Enter session numbers (based on those included in SPM.mat)'}';
