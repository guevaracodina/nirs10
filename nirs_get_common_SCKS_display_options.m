function DO = nirs_get_common_SCKS_display_options(job)
%Display options
DO.only_display = job.SCKS_display_options.only_display;
DO.show_mse = job.SCKS_display_options.show_mse;
DO.plot_algebraic_CMRO2 = job.SCKS_display_options.plot_algebraic_CMRO2;
DO.save_figures = job.SCKS_display_options.save_figures;
DO.generate_figures = job.SCKS_display_options.generate_figures;
DO.show_normalized_parameters = job.SCKS_display_options.show_normalized_parameters;
