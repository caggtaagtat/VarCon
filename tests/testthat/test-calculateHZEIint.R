test_that("calculateHZEIperNT works for example sequence", {
  x <- calculateHZEIperNT('ATACCAGCCAGCTATTACATTT')
  expect_equal(sum(x[,3]), 24.59)
})


test_that("calculateHZEIperNT works for example sequence", {
  x <- calculateHZEIperNT('ATACCAGCCAGCTATTACATTT')
  expect_is(x, "data.frame")
})

