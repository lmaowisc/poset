

## test if the wprod function works
test_that("wprod works", {
  expect_equal(wprod(1:3,rep(1,3)), 1)
  expect_equal(wprod(1:3,3:1), 0)
  expect_equal(wprod(1:3,3:5), -1)
})

