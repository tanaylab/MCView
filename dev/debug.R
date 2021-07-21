devtools::load_all()
init_config(project = project)
load_all_data(cache_dir = project_cache_dir(project))
init_defs()

future::plan(future::multicore)

options(golem.app.prod = FALSE)

shiny::devmode()