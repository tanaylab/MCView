# Run the MCView Application

Load project cache and run the MCView application.

## Usage

``` r
run_app(
  project,
  port = NULL,
  host = NULL,
  launch.browser = FALSE,
  profile = FALSE,
  ...
)
```

## Arguments

- project:

  path of the project to run.

- port:

  app port

- host:

  app host

- launch.browser:

  launch web browser after app start

- profile:

  enable profiling for debugging

- ...:

  Arguments passed on to
  [`shiny::shinyApp`](https://rdrr.io/pkg/shiny/man/shinyApp.html)

  `ui`

  :   The UI definition of the app (for example, a call to `fluidPage()`
      with nested controls).

      If bookmarking is enabled (see `enableBookmarking`), this must be
      a single argument function that returns the UI definition.

  `server`

  :   A function with three parameters: `input`, `output`, and
      `session`. The function is called once for each session ensuring
      that each app is independent.

  `onStart`

  :   A function that will be called before the app is actually run.
      This is only needed for `shinyAppObj`, since in the `shinyAppDir`
      case, a `global.R` file can be used for this purpose.

  `options`

  :   Named options that should be passed to the `runApp` call (these
      can be any of the following: "port", "launch.browser", "host",
      "quiet", "display.mode" and "test.mode"). You can also specify
      `width` and `height` parameters which provide a hint to the
      embedding environment about the ideal height/width for the app.

  `uiPattern`

  :   A regular expression that will be applied to each `GET` request to
      determine whether the `ui` should be used to handle the request.
      Note that the entire request path must match the regular
      expression in order for the match to be considered successful.

  `enableBookmarking`

  :   Can be one of `"url"`, `"server"`, or `"disable"`. The default
      value, `NULL`, will respect the setting from any previous calls to
      [`enableBookmarking()`](https://rdrr.io/pkg/shiny/man/enableBookmarking.html).
      See
      [`enableBookmarking()`](https://rdrr.io/pkg/shiny/man/enableBookmarking.html)
      for more information on bookmarking your app.

## Details

The app caches the plots and tables it computes (the shiny cache), by
default in memory, so a restarted app computes them again. When the
project's config sets `shiny_cache_dir` to a path, or when it doesn't
set it and the environment variable `MCVIEW_SHINY_CACHE_DIR` is set
(e.g. for all the apps of a shiny server), the cache is also kept on
disk, in a directory under that path that the next processes of the app
reuse until the project changes. Directories under that path that
weren't used for 30 days are removed. Apps that set
`MCVIEW_CLUSTERING_CACHE_DIR` share the clusterings of their heatmaps in
that directory (see
[`update_default_clusterings`](update_default_clusterings.md)).

## Examples

``` r
if (FALSE) { # \dontrun{
MCView::run_app("PBMC")
MCView::run_app(project = "PBMC", port = 5555, host = "127.0.0.1")
} # }
```
