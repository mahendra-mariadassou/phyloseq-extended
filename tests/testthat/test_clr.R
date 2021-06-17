test_that("clr throws an error for non positive data", {
  expect_error(clr(c(0, 1, 2)))
})

test_that("clr returns expected value for toy dataset", {
  expect_equal(clr(c(1, 2, 4)),
               c(-log(2), 0, log(2)))
})
