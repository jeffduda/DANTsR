context("antsrMesh")

test_that("isInImage returns TRUE", {
  x = makeImage(c(5,5))
  expect_true( isInImage(x, c(1,1,1)) )
})

test_that("isInImage returns FALSE", {
  x = makeImage(c(5,5))
  expect_false( isInImage(x, c(-1,1,1)) )
})

test_that("pointIsInImage returns TRUE", {
  x = makeImage(c(5,5))
  expect_true( pointIsInImage(x, c(1,1,1)) )
})

test_that("pointIsInImage returns FALSE", {
  x = makeImage(c(5,5))
  expect_false( pointIsInImage(x, c(-1,1,1)) )
})

test_that("indexIsInImage returns TRUE", {
  x = makeImage(c(5,5))
  expect_true( indexIsInImage(x, c(1,1,1)) )
})

test_that("indexIsInImage returns FALSE", {
  x = makeImage(c(5,5))
  expect_false( indexIsInImage(x, c(-1,1,1)) )
})
