devtools::load_all()
project <- get_default_project()
init_config(project = project)
load_all_data(cache_dir = project_cache_dir(project))
init_defs()

future::plan(future::multicore)