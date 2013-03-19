function two_sample_t_test = nirs_dfg_two_sample_t_test

group1         = cfg_entry;
group1.tag    = 'group1';
group1.name     = 'Subjects for group 1';
group1.strtype = 'r';
group1.num     = [1 Inf];
group1.help    = {'Enter the subject numbers for group 1.'
    'This is according to the order of the NIRS.mat structures passed from module to module.'}';

group2         = cfg_entry;
group2.tag    = 'group2';
group2.name     = 'Subjects for group 2';
group2.strtype = 'r';
group2.num     = [1 Inf];
group2.help    = {'Enter the subject numbers for group 2.'};

eq_variance = cfg_menu;
eq_variance.tag  = 'eq_variance';
eq_variance.name = 'Assume equal variance';
eq_variance.labels = {'Yes','No'};
eq_variance.values = {1,0};
eq_variance.val = {1};
eq_variance.help = {'Assume that the variance of group 2 is the same as the variance of group 1.'};

two_samples         = cfg_branch;
two_samples.tag     = 'two_samples';
two_samples.name    = 'Two samples'; 
two_samples.val     = {group1 group2 eq_variance};
two_samples.help    = {'Two samples'}';

one_sample         = cfg_branch;
one_sample.tag     = 'one_sample';
one_sample.name    = 'One sample'; 
one_sample.val     = {};
one_sample.help    = {'One sample.'}';

two_sample_t_test           = cfg_choice;
two_sample_t_test.name      = 'Choose contrast choice method';
two_sample_t_test.tag       = 'two_sample_t_test';
two_sample_t_test.values    = {one_sample two_samples}; 
two_sample_t_test.val       = {one_sample}; 
two_sample_t_test.help      = {'Choose to do 1 sample or 2 samples t test'}'; 
