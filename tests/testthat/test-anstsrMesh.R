context("antsrMesh")

test_that("mesh can be created with new()", {
  x = new("antsrMesh")
  expect_true( class(x)=="antsrMesh")
})
