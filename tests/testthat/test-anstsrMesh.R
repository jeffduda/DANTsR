context("antsrMesh")

test_that("mesh can be created with new()", {
  x = new("antsrMesh")
  expect_true( class(x)=="antsrMesh")
})

test_that("mesh can be created with new()", {
  x = new("antsrMesh", reserve=10)
  expect_true( ( class(x)=="antsrMesh" ) &
               ( antsrMeshGetNumberOfPoints(x)==10 ) )
})


# Test all parameter combos for the C++ templated parameters
for ( d in c(2,3,4) ) {
  for ( p in c("float", "double") ) {


    test_that("mesh can be created with new()", {
      x = new("antsrMesh", dimension=d, precision=p)
      expect_true( ( class(x)=="antsrMesh" ) &
                   ( x@dimension == d ) &
                   ( x@precision == p ) )
    })
  }
}
