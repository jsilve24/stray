D <- 6
P <- 4
A.prop <- array(abs(rnorm(10*6*5)), dim=c(10, 6, 5))
A.prop[1:D,,] <- driver::miniclo_array(A.prop[1:D,,], parts=1)

test_that("orthus data transform correctness", {
  # TEST ilr and ilrinv
  A.ilr <- oilr(A.prop, D)
  expect_equal(A.ilr[1:(D-1),,], ilr_array(A.prop[1:D,,],parts=1))
  expect_equal(A.prop[(D+1):(D+P),,], A.ilr[(D):(D-1+P),,])
  expect_equal(oilrInv(A.ilr,D-1), A.prop)
  rm(A.ilr)
  
  # TEST alr and alrinv
  A.alr <- oalr(A.prop, D)
  expect_equal(A.alr[1:(D-1),,], alr_array(A.prop[1:D,,],parts=1))
  expect_equal(A.prop[(D+1):(D+P),,], A.alr[(D):(D-1+P),,])
  expect_equal(oalrInv(A.alr,D-1), A.prop)
  rm(A.alr)
  
  # TEST clr and clrinv
  A.clr <- oclr(A.prop, D)
  expect_equal(A.clr[1:(D),,], clr_array(A.prop[1:D,,],parts=1))
  expect_equal(A.prop[(D+1):(D+P),,], A.clr[(D+1):(D+P),,])
  expect_equal(oclrInv(A.clr,D), A.prop)
  rm(A.clr)
})


Sigma.ilr <- rWishart(5, 20, diag(D-1+P))

test_that("orthus covariance transform correctness", {
  # test oilrvar2clrvar
  V <- create_default_ilr_base(D)
  Sigma.clr <- oilrvar2clrvar(Sigma.ilr, D-1, V)
  expect_equal(Sigma.clr[1:D,1:D,1], ilrvar2clrvar(Sigma.ilr[1:(D-1),1:(D-1),1],V))
  foo <- ilrInv_array(Sigma.ilr[1:(D-1), D:(D-1+P),],coords=1)
  foo <- clr_array(foo, parts=1)
  expect_equal(Sigma.clr[1:D,(D+1):(D+P),], foo)
  expect_equal(Sigma.clr[(D+1):(D+P),1:D,], aperm(foo, c(2,1,3)))
  expect_equal(Sigma.clr[(D+1):(D+P),(D+1):(D+P),], Sigma.ilr[(D):(D+P-1),(D):(D+P-1),])
  
  # test oclrvar2ilrvar
  expect_equal(Sigma.ilr, oclrvar2ilrvar(Sigma.clr, D, V))
  
  # test oalrvar2clrvar
  Sigma.alr <- Sigma.ilr; 
  Sigma.clr <- oalrvar2clrvar(Sigma.alr, D-1, D)
  expect_equal(Sigma.clr[1:D,1:D,1], alrvar2clrvar(Sigma.alr[1:(D-1),1:(D-1),1],D))
  foo <- alrInv_array(Sigma.alr[1:(D-1), D:(D-1+P),],coords=1)
  foo <- clr_array(foo, parts=1)
  expect_equal(Sigma.clr[1:D,(D+1):(D+P),], foo)
  expect_equal(Sigma.clr[(D+1):(D+P),1:D,], aperm(foo, c(2,1,3)))
  expect_equal(Sigma.clr[(D+1):(D+P),(D+1):(D+P),], Sigma.alr[(D):(D+P-1),(D):(D+P-1),])
  
  # test oclrvar2alrvar
  expect_equal(Sigma.alr, oclrvar2alrvar(Sigma.clr, D, D))
  
  # test oilrvar2ilrvar
  Sigma.ilr <- Sigma.alr
  expect_equal(Sigma.ilr, oilrvar2ilrvar(Sigma.ilr, D-1, V, V))
  
  # Others not currently tested as they are just based on the above transforms
})


test_that("orthusfit transforms don't give error", {
  sim <- orthus_sim()
  fit <- orthus(sim$Y, sim$Z, sim$X)
  fit <- to_proportions(fit)
  fit <- to_alr(fit, 4)
  fit <- to_ilr(fit)
  fit <- to_clr(fit)
  expect_true(TRUE)
})


