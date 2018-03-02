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

test_that("mesh can be created with antsrMeshCreate()", {
  x = antsrMeshCreate()
  expect_true( class(x)=="antsrMesh")
})

test_that("mesh can be created with antsrMeshCreate()", {
  x = antsrMeshCreate(reserve=10)
  expect_true( ( class(x)=="antsrMesh" ) &
               ( antsrMeshGetNumberOfPoints(x)==10 ) )
})

# Test all parameter combos for the C++ templated parameters
for ( d in c(2,3,4) ) {
  for ( p in c("float", "double") ) {

    test_that("mesh can be created with antsrMeshCreate()", {
      x = antsrMeshCreate(dimension=d, precision=p)
      expect_true( ( class(x)=="antsrMesh" ) &
                   ( x@dimension == d ) &
                   ( x@precision == p ) )
    })
  }
}

test_that( "point can be added to a mesh with no index", {
  x = antsrMeshCreate()
  antsrMeshAddPoint(x, c(1,2,3) )
  antsrMeshAddPoint(x, c(1,2,3) )
  expect_true( antsrMeshGetNumberOfPoints(x) == 2)
})

test_that( "point can be retrieved", {
  x = antsrMeshCreate()
  antsrMeshAddPoint(x, c(1,2,3) )
  expect_true( sum(antsrMeshGetPoint(x,0) == c(1,2,3))==3 )
})

test_that( "point can be added to a mesh with an index", {
  x = antsrMeshCreate()
  antsrMeshAddPoint(x, c(1,2,3), 5 )
  expect_true( antsrMeshGetNumberOfPoints(x) == 6)
})

test_that( "points can be retrieved as matrix", {
  x = antsrMeshCreate()
  for (i in 1:4) {
    antsrMeshAddPoint(x, i*c(1,2,3), i-1 )
  }
  pts = antsrMeshGetPoints(x)
  expect_true( ( class(pts) == "matrix" ) &
               ( sum( dim(pts) == c(4,3) ) == 2 ) &
               ( sum(pts) == 60 ) )
})

#test_that("matrix of points can be passed to mesh", {
#  x = antsrMeshCreate()
#  m = matrix(1:12, 4, 3)
#  antsrMeshSetPoints(x, m)
#}
