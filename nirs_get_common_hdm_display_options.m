function DO = nirs_get_common_hdm_display_options(job)
%Display options
DO.only_display = job.hdm_display_options.only_display;
DO.show_mse = job.hdm_display_options.show_mse;
DO.plot_algebraic_CMRO2 = job.hdm_display_options.plot_algebraic_CMRO2;
DO.save_figures = job.hdm_display_options.save_figures;
DO.generate_figures = job.hdm_display_options.generate_figures;
DO.show_normalized_parameters = job.hdm_display_options.show_normalized_parameters;
