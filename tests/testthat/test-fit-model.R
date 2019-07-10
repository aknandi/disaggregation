
context("Fitting model")

test_that("fit_model produces errors whe expected", {
  
  test_data <- list(polygon_data = data.frame(),
                   covariate_data = data.frame(),
                   coords = matrix(),
                   startendindex = matrix(),
                   mesh = NULL)
  
  expect_error(fit_model(list()))
  
})