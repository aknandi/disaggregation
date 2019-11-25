library(testthat)
library(disaggregation)

not_cran <- function() identical(Sys.getenv("NOT_CRAN"), "true")

if(not_cran()) {
  test_check("disaggregation")
}
