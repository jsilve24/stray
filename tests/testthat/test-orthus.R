test_that("orthus sim and wrapper run without error", {
  sim <- orthus_sim()
  fit <- orthus(sim$Y, sim$Z, sim$X)
  expect_true(TRUE)
})
