devtools::load_all()
project <- get_default_project()
init_config(project = project)
load_all_data(data_dir = project_data_dir(project))
init_defs()

future::plan(future::multicore)