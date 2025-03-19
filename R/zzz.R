globalVariables(unique(c(
    # add_gene_modules:
    "gene", "module",
    # add_genes_to_gene_module:
    "gene",
    # add_markers_colorbars:
    "metacell",
    # app_server:
    "cell_type", "color", "config", "module", "project", "tab_defs",
    # app_server : dataset:
    "project",
    # app_ui:
    "project", "tab_defs",
    # atlas_gene_gene:
    "default_gene1", "default_gene2", "egc_epsilon",
    # calc_ct_ct_gene_df:
    "egc_epsilon",
    # calc_gene_modules:
    "clust", "gene", "module",
    # calc_gg_mc_top_cor:
    "col1", "col2", "val",
    # calc_mc_mc_gene_df:
    "egc_epsilon",
    # calc_obs_exp_mc_df:
    "egc_epsilon",
    # calc_obs_exp_type_df:
    "egc_epsilon",
    # calc_samp_mc_count:
    "mc_data", "metacell", "samp_id",
    # calc_samp_mc_frac:
    "mc_data",
    # calc_samp_metadata:
    "mc_data", "samp_id",
    # calc_samp_samp_gene_df:
    "egc_epsilon",
    # calc_samples_list:
    "mc_data",
    # cell_metadata_to_metacell:
    "cat_var_n", "cat_varp", "cell_id", "metacell",
    # cell_metadata_to_metacell_from_h5ad:
    "cell_id", "metacell_name", "outlier",
    # cell_type_gene_boxplot:
    "Cell type", "cell_type", "egc_epsilon", "expr_breaks", "metacell",
    # cell_type_metadata_boxplot:
    "Cell type", "cell_type", "metacell",
    # cell_type_metadata_confusion:
    "# of metacells", "% of cell type metacells", "Cell type", "cell_type", "color", "metacell", "n_tot", "n_tot_md", "p_cell_type", "p_md", "total # of cell type metacells",
    # choose_markers:
    "fp", "gene", "metacell",
    # clipboard_reactives:
    "cell_type", "metacell",
    # clipboard_reactives : <anonymous>:
    "cell_type", "metacell",
    # colored_metacell_selector:
    "config",
    # compute_umap:
    "idx", "k", "mc1", "mc2", "weight",
    # config_shiny_cache:
    "config", "project",
    # download_modal_reactives : <anonymous>:
    "project",
    # dt_selector_column:
    "cell_type",
    # filter_mat_by_cell_types:
    "cell_type", "metacell",
    # fitted_genes_per_cell_type_plot:
    "cell_type", "color", "type",
    # fitted_genes_per_cell_type_table:
    "common", "fitted_gene_of", "gene", "type",
    # gene_correction_factor_scatter_plot:
    "Correction factor", "correction_factor", "gene", "Gene", "Max expression", "max_expr", "type", "Type",
    # gene_correction_factor_table:
    "correction_factor", "gene", "Max expression", "max_expr", "type",
    # gene_inner_fold_scatter_plot:
    "# of metacells", "gene", "Gene", "Max expression", "max_expr", "significant_inner_folds_count", "type", "Type",
    # gene_inner_fold_table:
    "# of metacells", "gene", "Max expression", "max_expr", "significant_inner_folds_count", "type",
    # gene_label:
    "gene", "module",
    # get_cell_type_colors:
    "cell_type", "color",
    # get_cell_type_data:
    "cell_type", "mc_data",
    # get_cell_types_mat:
    "cell_type", "metacell",
    # get_marker_genes:
    "config",
    # get_markers:
    "cache_dir",
    # get_mc_color_key:
    "cell_type", "mc_col",
    # get_mc_config:
    "config",
    # get_mc_data:
    "mc_data",
    # get_mc_gene_modules_egc:
    "module",
    # get_metacell_types_data:
    "cell_type", "color", "mc_data", "metacell",
    # get_metadata:
    "mc_data",
    # get_module_genes:
    "gene",
    # get_samples_egc:
    "cell_type", "metacell",
    # get_samples_mat:
    "cell_type", "metacell",
    # get_top_cor_gene:
    "egc_epsilon", "gene1", "gene2", "label",
    # golem_add_external_resources:
    "config",
    # group_selectors:
    "cell_type", "metacell",
    # group_selectors_mod_query:
    "cell_type", "metacell",
    # has_atlas:
    "mc_data",
    # heatmap_matrix_reactives:
    "cell_type", "config", "metacell",
    # heatmap_reactives : <anonymous>:
    "cell_type", "color", "type",
    # heatmap_sidebar:
    "config",
    # help_reactives:
    "config", "help_config",
    # import_atlas:
    "atlas", "cell_type", "color", "dtype", "gene", "metacell", "projected_type", "similar", "top1_gene", "top1_lfp", "top2_gene", "top2_lfp", "type", "value", "weight",
    # import_cell_metadata:
    "cell_id", "metacell", "samp_id",
    # import_dataset:
    "cell_type", "cluster", "color", "correction_factor", "fp", "gene", "gene1", "gene2", "grouped", "metacell", "significant_inner_folds_count", "total_umis",
    # import_dataset_metacell1:
    "cell_type", "color", "metacell", "time_bin",
    # init_config:
    "config",
    # init_selected_genes:
    "config", "project",
    # init_tab_defs:
    "config", "tab_defs",
    # init_temp_scdb:
    ".",
    # initial_proj_point_size:
    "config",
    # initial_proj_stroke:
    "config",
    # initial_scatters_point_size:
    "config",
    # initial_scatters_stroke:
    "config",
    # layout_and_graph_to_mc2d:
    "from", "metacell", "to",
    # load_all_data:
    "project",
    # load_all_mc_data:
    ".", "mc_data",
    # load_all_mc_data_atlas:
    ".", "mc_data",
    # load_default_2d_projection:
    "i", "j", "x", "y",
    # load_outliers:
    "fp",
    # mc2d_plot_gene_ggp:
    "cell_type", "egc_epsilon", "enrich", "expr_clipped", "mc_age", "metacell",
    # mc2d_plot_metadata_ggp:
    "Cell type", "cell_type", "mc_age", "mc_col", "metacell", "top1_gene", "top1_lfp", "top2_gene", "top2_lfp",
    # mc2d_to_graph_df:
    "d", "dx", "dy", "from", "mc1", "mc2", "metacell", "range_x", "range_y", "to", "weight", "x", "x_mc1", "x_mc2", "y", "y_mc1", "y_mc2",
    # mc_mc_gene_scatter_df_reactive:
    "cell_type", "metacell",
    # mctnetwork_g_t_types:
    "flow", "g_src", "g_targ", "mc_t1", "mc_t2", "src_g", "targ_g", "time1", "time2", "type1", "type2",
    # metacell_colors_reactive:
    "mc_col", "metacell",
    # metacell_from_coords_proj:
    "diff1", "diff2",
    # metacell_selector:
    "config",
    # metacell_selectors:
    "config",
    # metacell_selectors_mod_query:
    "cell_type", "config",
    # metacells_to_types:
    "cell_type", "metacell",
    # min_edge_length:
    "config",
    # mod_about_server : <anonymous>:
    "cell_type", "color", "default_gene1", "default_gene2", "x", "x_mc1", "x_mc2", "y", "y_mc1", "y_mc2",
    # mod_about_ui:
    "about_file", "config",
    # mod_annotate_server : <anonymous>:
    "color", "metacell", "project",
    # mod_annotate_server : <anonymous> : <anonymous>:
    "color", "metacell", "top1_gene", "top1_lfp", "top2_gene", "top2_lfp",
    # mod_atlas_server : <anonymous>:
    "project",
    # mod_cell_type_server : <anonymous>:
    "default_gene1", "egc_epsilon",
    # mod_flow_server : <anonymous>:
    "cell_type", "color", "default_gene1", "metacell",
    # mod_gene_module_controllers:
    "gene", "module",
    # mod_gene_module_controllers : <anonymous>:
    "gene", "module",
    # mod_manifold_server : <anonymous>:
    "mc1", "mc2", "project",
    # mod_manifold_server : <anonymous> : <anonymous>:
    "mc1", "mc2",
    # mod_outliers_server : <anonymous>:
    "cell_id", "metacell", "most_similar",
    # mod_projection_qc_server : <anonymous>:
    "fitted_gene_any", "gene", "marker_gene",
    # mod_query_server : <anonymous>:
    "cell_type", "default_gene1", "fitted_gene", "fitted_gene_any", "marker_gene", "metacell",
    # mod_samples_server : <anonymous>:
    "default_gene1", "default_gene2", "samp_id",
    # modify_gene_names:
    "alt_name", "gene_name", "new_name",
    # observe_mc_click_event:
    "cell_type", "metacell",
    # observer_mc_select_event:
    "cell_type", "metacell",
    # order_mc_by_most_var_genes:
    ".", "cell_type", "ct_ord", "glob_ord", "metacell", "orig_ord", "rand",
    # parse_cell_type_colors:
    "cell_type", "cluster", "color",
    # parse_metacell_types:
    "age", "cluster", "metacell",
    # parse_metadata:
    "metacell",
    # plot_gene_time_over_mc:
    "Age", "Cell type", "cell_type", "egc_epsilon", "expr_breaks", "mc_age", "metacell", "Metacell",
    # plot_gene_trajectory:
    "egc_epsilon", "expr_breaks", "gene", "Gene", "time_bin", "time_desc",
    # plot_gg_over_mc:
    "Cell type", "cell_type", "egc_epsilon", "expr_breaks", "mc_age", "metacell", "Metacell",
    # plot_markers_mat:
    "cell_type", "color",
    # plot_mc_mc_gene:
    "expr_breaks", "gene", "Gene", "pval",
    # plot_mc_scatter:
    "cell_type", "color", "color_values", "egc_epsilon", "expr_breaks", "metacell", "Metacell",
    # plot_mc_time_dist:
    "Metacell age (E[t])", "n_cells", "time_bin", "time_desc",
    # plot_obs_proj_scatter:
    "cell_type", "color", "color_values", "egc_epsilon", "expr_breaks", "metacell", "Metacell", "query", "weight",
    # plot_propagation_net_metacell:
    "cell_type", "cell_type1", "edge_lwd", "edge_lwd_fg", "edge_type", "flow", "mc_col", "mc1", "mc2", "time1", "time2", "type1", "type2", "x1", "x2", "y1", "y2",
    # plot_sample_scatter:
    "cell_type", "color", "color_values", "egc_epsilon", "expr_breaks", "frac", "metacell", "samp_id", "Sample",
    # plot_sample_scatter : get_var_md:
    "samp_id",
    # plot_sample_stacked_types:
    "# of cells", "Cell type", "cell_type", "color", "Fraction", "metacell", "samp_id", "Sample",
    # plot_type_predictions_bar:
    "cell_type", "color", "fraction", "metacell", "type",
    # plot_vein:
    ".",
    # projection_selectors:
    "default_gene1",
    # qc_ecdf:
    "x", "y",
    # read_metacell_graphs : <anonymous>:
    "from", "to", "weight",
    # render_2d_plotly:
    "cell_type", "metacell", "query", "weight",
    # render_2d_plotly : plot_2d_atlas_proj:
    "cell_type", "metacell", "query", "weight", "Weight",
    # render_mc_mc_gene_diff_table:
    "gene", "Gene", "pval", "type",
    # render_mc_mc_gene_plotly:
    "cell_type", "fitted_gene",
    # scatter_box_outputs:
    "default_gene1", "default_gene2", "egc_epsilon", "metacell",
    # select_top_fold_genes_per_metacell:
    "fp", "gene", "metacell",
    # top_correlated_selector:
    "project",
    # top_correlated_selector_multiple_genes:
    "project",
    # update_metacell_types:
    "cell_type", "metacell",
    # update_metadata:
    "metacell",
    # verify_gene_modules:
    "gene",
    # zero_fold_gene_plot:
    "avg", "Expected", "Expression", "FC", "gene", "Gene", "metacell", "Metacell", "obs", "Observed", "type", "Type", "zero_fold",
    # zero_fold_table:
    "avg", "Expected", "Expression", "FC", "gene", "metacell", "obs", "type", "zero_fold",
    ":=", "sec_ord", "metacell_names", "scdb_mctnetwork", "scdb_mc", "mctnetwork_get_flow_mat"
)))
