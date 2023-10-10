
context("Build mesh")

test_that("build_mesh behaves as expected", {

  skip_on_cran()

  my_mesh <- build_mesh(spdf)

  expect_error(build_mesh(response_df))
  expect_error(build_mesh(spdf, mesh.args = c(4, 8)))
  expect_is(my_mesh, 'inla.mesh')
  expect_is(build_mesh(spdf, mesh.args = list(max.edge = c(50, 100))), 'inla.mesh')

})
