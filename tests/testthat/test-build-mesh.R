
context("Build mesh")

test_that("build_mesh behaves as expected", {

  my_mesh <- build_mesh(spdf)

  expect_error(build_mesh(response_df))
  expect_error(build_mesh(spdf, mesh_args = c(4, 8)))
  expect_message(build_mesh(spdf, mesh_args = list(cut = 0.1)))
  expect_is(my_mesh, 'inla.mesh')
  expect_is(build_mesh(spdf, mesh_args = list(max.edge = c(50, 100))), 'inla.mesh')

})
