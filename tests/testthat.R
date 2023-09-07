library(testthat)
library(MCView)

if (identical(Sys.getenv("NOT_CRAN"), "true")) {
    test_check("MCView")
}
