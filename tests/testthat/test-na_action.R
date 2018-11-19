context("NA actions")
set.seed(12345)
Y <- matrix(rnorm(100*20), nrow = 100)
Y_row <- Y
Y_row[5,] <- NA_real_

Y_random <- Y
cells_samp <- cbind(sample(n, size = p, replace = TRUE),
                    sample(p, size = p, replace = TRUE))
for (i in seq_len(nrow(cells_samp))) {
  Y_random[cells_samp[i,1],
           cells_samp[i,2]] <- NA_real_
}
X <- rnorm(100)

test_that("NA action is fail", {
  pcev_out <- try(computePCEV(Y, X,
                              na_action = "fail"),
                  silent = TRUE)
  pcev_out2 <- try(computePCEV(Y_row, X,
                               na_action = "fail"),
                   silent = TRUE)
  pcev_out3 <- try(computePCEV(Y_random, X,
                               na_action = "fail"),
                   silent = TRUE)
  
  expect_false(inherits(pcev_out, "try-error"))
  expect_true(inherits(pcev_out2, "try-error"))
  expect_true(inherits(pcev_out3, "try-error"))
})

test_that("NA action is omit", {
  pcev_out <- try(computePCEV(Y, X,
                              na_action = "omit"),
                  silent = TRUE)
  pcev_out2 <- try(computePCEV(Y_row, X,
                               na_action = "omit"),
                   silent = TRUE)
  pcev_out3 <- try(computePCEV(Y_random, X,
                               na_action = "omit"),
                   silent = TRUE)
  
  expect_false(inherits(pcev_out, "try-error"))
  expect_false(inherits(pcev_out2, "try-error"))
  expect_false(inherits(pcev_out3, "try-error"))
})

test_that("NA action is column", {
  pcev_out <- try(computePCEV(Y, X,
                              na_action = "column"),
                  silent = TRUE)
  pcev_out2 <- try(computePCEV(Y_row, X,
                               na_action = "column"),
                   silent = TRUE)
  pcev_out3 <- try(computePCEV(Y_random, X,
                               na_action = "column"),
                   silent = TRUE)
  
  expect_false(inherits(pcev_out, "try-error"))
  expect_false(inherits(pcev_out2, "try-error"))
  expect_false(inherits(pcev_out3, "try-error"))
})
