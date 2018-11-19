context("Base tests")

Y <- matrix(rnorm(100*20), nrow = 100)
X <- rnorm(100)

test_that("No errors with classical", {
  pcev_out <- try(computePCEV(Y, X),
                  silent = TRUE)
  pcev_out2 <- try(computePCEV(Y, X, shrink = TRUE),
                   silent = TRUE)
  
  expect_false(inherits(pcev_out, "try-error"))
  expect_false(inherits(pcev_out2, "try-error"))
})

test_that("No errors with singular", {
  pcev_out <- try(computePCEV(Y, X, estimation = "singular",
                              nperm = 5),
                  silent = TRUE)
  pcev_out2 <- try(computePCEV(Y, X, estimation = "singular",
                               nperm = 5, shrink = TRUE),
                   silent = TRUE)
  
  expect_false(inherits(pcev_out, "try-error"))
  expect_false(inherits(pcev_out2, "try-error"))
})

test_that("No errors with block", {
  pcev_out <- try(computePCEV(Y, X, estimation = "block",
                              inference = "permutation",
                              nperm = 5),
                  silent = TRUE)
  pcev_out2 <- try(computePCEV(Y, X, estimation = "block",
                               shrink = TRUE,
                               inference = "permutation",
                               nperm = 5),
                   silent = TRUE)
  
  expect_false(inherits(pcev_out, "try-error"))
  expect_false(inherits(pcev_out2, "try-error"))
})
