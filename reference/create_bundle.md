# Generate a 'deployment ready' bundle of the a project app

Generate a 'deployment ready' bundle of the a project app

## Usage

``` r
create_bundle(
  project,
  path = getwd(),
  name = "MCView_bundle",
  overwrite = FALSE,
  self_contained = FALSE,
  branch = "latest_release",
  restart = overwrite,
  permissions = NULL,
  light_version = FALSE,
  excluded_tabs = c("Gene modules", "Annotate", "Inner-fold", "Stdev-fold",
    "Gene correlation"),
  shiny_cache_dir = NULL,
  shiny_cache_max_size = NULL,
  ...
)
```

## Arguments

- project:

  path to the project directory

- path:

  path in which to create the bundle.

- name:

  name of the folder in which to create the bundle. The bundle would be
  created at `path`/`name`

- overwrite:

  overwrite bundle if already exists

- self_contained:

  include the source code of `MCView` in the bundle and use it to run
  the app. Use this in order to ensure that the package would always run
  the same way, regardless of MCView changes. When this option is FALSE,
  the version of `MCView` which is installed on the server would be
  loaded, which can be occasionally be different than the one used when
  creating the app. By default, the code uses the latest `MCView`
  release would be used, see `branch` for other options.

- branch:

  name of the `MCView` branch to include when `self_contained=TRUE`. By
  default, the latest release would be used. You can set this parameter
  to NULL in order to include the current development version ('master'
  branch), or set it to any other branch in the 'tanaylab/MCView' github
  repository.

- restart:

  add a file named 'restart.txt' to the bundle. This would force
  shiny-server to restart the app when updated.

- permissions:

  change the file permissions of the bundle after creation, e.g. "777".
  When NULL - permissions would not be changed.

- light_version:

  create a light version of the bundle, which would not include features
  that require heavy computation (e.g. changing Marker genes, Gene
  modules etc.)

- excluded_tabs:

  a character vector of tabs to exclude from the light version of the
  bundle.

- shiny_cache_dir:

  a path to a directory in which to store shiny cache, can be relative
  to the bundle, e.g. "./shiny_cache". If set to TRUE, a temporary
  directory would be set. If NULL - shiny would cache objects in memory.

- shiny_cache_max_size:

  maximum size of the shiny cache in bytes. Default is 200e6.

- ...:

  Arguments passed on to
  [`gert::git_clone`](https://docs.ropensci.org/gert/reference/git_fetch.html)

  `password`

  :   a string or a callback function to get passwords for
      authentication or password protected ssh keys. Defaults to
      [askpass](https://r-lib.r-universe.dev/askpass/reference/askpass.html)
      which checks `getOption('askpass')`.

  `ssh_key`

  :   path or object containing your ssh private key. By default we look
      for keys in `ssh-agent` and
      [credentials::ssh_key_info](https://docs.ropensci.org/credentials/reference/ssh_credentials.html).

  `verbose`

  :   display some progress info while downloading

  `mirror`

  :   use the `--mirror` flag

  `url`

  :   remote url. Typically starts with `https://github.com/` for public
      repositories, and `https://yourname@github.com/` or
      `git@github.com/` for private repos. You will be prompted for a
      password or pat when needed.

  `bare`

  :   use the `--bare` flag

## Details

Create a minimal shiny app in `path`/`name` directory which would
contain:

- app.R file.

- project config and cache.

The bundle can then be deployed in shiny-server, shinyapps.io or any
other environment that supports serving shiny apps.

Note: when deploying to these services - make sure you have the MCView
package installed.

## Examples

``` r
if (FALSE) { # \dontrun{
MCView::create_bundle(project = "PBMC", path = getwd(), name = "PBMC")

# latest release
MCView::create_bundle(project = "PBMC", path = getwd(), name = "PBMC", self_contained = TRUE)

# development version
MCView::create_bundle(
    project = "PBMC",
    path = getwd(),
    name = "PBMC",
    self_contained = TRUE,
    branch = NULL
)

# specific branch
MCView::create_bundle(
    project = "PBMC",
    path = getwd(),
    name = "PBMC",
    self_contained = TRUE,
    branch = "feat@atlas-projection"
)
} # }
```
