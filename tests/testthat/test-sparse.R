library(Matrix)
library(Rvolesti)

test_that("Sparse Hpolytope works", {

  A <- rsparsematrix(10, 5, density = 0.2)
  b <- rep(1, 10)

  P <- Hpolytope(A, b)

  expect_true(inherits(P, "Hpolytope_sparse"))
})