context("test-random_pibble_init.R")

test_that("random_pibble_init works", {
  Y <- matrix(sample(1:100, 100), 10, 10)
  foo <- random_pibble_init(Y)
  expect_equal(dim(foo), c(9,10))
})


test_that("check_dims correct", {
  y <- c(1,2,3)
  expect_error(expect_error(check_dims(y, 3, "y"))) # expect no error!
  expect_error(check_dims(y, c(3,1), "y"))
  
  y <- matrix(c(1,3,4,5), 2, 2)
  expect_error(expect_error(check_dims(y, c(2,2), "y"))) #expect no error!
  expect_error(check_dims(y, c(2), "y"))
})

# test_that("name correct on unnamed imput", {
#   sim <- pibble_sim()
#   sim$Y <- unname(sim$Y)
#   sim$X <- unname(sim$X)
#   attach(sim)
#   fit <- pibble(Y, X)
#   
#   
#   
#   # When not all parameters are present
#   fit$Eta <- NULL
#   name(fit)
#   
#   detach(sim)
#   expect_true(TRUE) # so that above does not give error
# })