test_that("clr throws an error for non positive data", {
  expect_error(clr(c(0, 1, 2)))
})

test_that("clr returns expected value for toy dataset", {
  expect_equal(clr(c(1, 2, 4)),
               c(-log(2), 0, log(2)))
})

test_that("mclr returns expected value for toy dataset", {
  shift <- log(2)/2 + 1
  expect_equal(mclr(matrix(c(0,1,2,1,1,0), nrow = 2, byrow = TRUE), c = 1),
               matrix(c(0, 1, log(2)/2 + shift, shift, shift, 0), nrow = 2, byrow = TRUE))
  shift <- log(2)/2 + 4
  expect_equal(mclr(matrix(c(0,1,2,1,1,0), nrow = 2, byrow = TRUE), c = 4),
               matrix(c(0, 4, log(2)/2 + shift, shift, shift, 0), nrow = 2, byrow = TRUE))
})
