
context("Match points to polygons")

test_that("Getting start and end index errors when it should", {

  covs <- data.frame(area_id = c(1, 1, 1, 2, 2, 3, 3, 3, 3), cov1 = c(3, 9, 5, 2, 3, 6, 7, 3, 5))
  response <- data.frame(area_id = c(1, 2, 3), response = c(4, 7, 2))

  expect_error(getStartendindex(covs, response, 'id'))
  expect_error(getStartendindex(response, covs, 'id'))
})

test_that("Getting start and end index returns the right object", {

  covs <- data.frame(area_id = c(1, 1, 1, 2, 2, 3, 3, 3, 3, 4), cov1 = c(3, 9, 5, 2, 3, 6, 7, 3, 5, 6))
  response <- data.frame(area_id = c(1, 2, 3), response = c(4, 7, 2))

  result <- getStartendindex(covs, response, 'area_id')

  expect_is(result, "matrix")
  expect_equal(nrow(result), nrow(response))
  expect_equal(ncol(result), 2)
  expect_equal(result, matrix(c(0, 3, 5, 2, 4, 8), nrow = 3, ncol = 2))
})
